# ========================
# PyRPOD: tests/mdao/mdao_unit_test_07.py
# ========================
# Unit tests for the ISS-panel study extensions:
#
#   * plume-model dispatch -- both collisionless Cai variants are selectable
#     by name, unknown names are rejected, and BOTH reduce to the same
#     LocalFieldState structure (the interface the shared Maxwellian
#     gas-surface formulas consume);
#   * panel-local pose generation -- a centered normal-incidence source
#     lands where it should with the expected DCM, a nonzero u offset
#     translates the source parallel to the panel WITHOUT tilting the plume
#     axis, and offsets default to zero;
#   * sweep enumeration -- the offset grid produces exactly
#     n_distances x n_u_offsets x n_v_offsets cases, in a deterministic
#     order, and incompatible pose definitions are refused;
#   * derived Knudsen metadata -- correct for both reference-length modes,
#     cleanly absent when unconfigured, and clearly rejected when invalid;
#   * panel-local coordinate transforms and load projections, including a
#     hand-computed moment and normal-force case and the sign conventions;
#   * the distribution-CSV schema.
#
# Nothing here touches DSMC: every model under test is collisionless, and
# the Knudsen number is asserted to be metadata that no solution reads.
#
# Run:  python -m pytest mdao/mdao_unit_test_07.py   (from tests/)

import copy
import csv
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pytest
import yaml

from pyrpod.mdao.firing_plan import (
    build_case_firings,
    pose_for,
    pose_for_sweep_pose,
    translated_pose_for,
)
from pyrpod.mdao.study_config import (
    KnudsenSpec,
    StudyConfig,
    StudyConfigError,
    SweepPose,
    SweepSpec,
    TargetSpec,
)
from pyrpod.mdao.study_results import CaseResult
from pyrpod.mdao.surface_distribution import (
    DISTRIBUTION_COLUMNS,
    distribution_rows,
    write_surface_distribution,
)
from pyrpod.mdao.surface_loads import (
    integrate_component_loads,
    panel_local_coordinates,
    project_to_panel_frame,
)
from pyrpod.plume.gas_kinetics_models import (
    DEFAULT_PLUME_MODEL,
    PLUME_MODELS,
    PlumeModelError,
    create_model,
    kinetics_key_for,
    local_field_state,
    maxwellian_surface_loads,
    model_name_for_kinetics,
    resolve_model_name,
)

_TESTS_DIR = Path(__file__).resolve().parents[1]
FLAT_PLATE_CASE = _TESTS_DIR.parent / 'case' / 'plume' / 'plume_flat_plate_sweep'
ISS_CASE = _TESTS_DIR.parent / 'case' / 'plume' / 'iss_panel_thesis'
BASELINE_YAML = FLAT_PLATE_CASE / 'study' / 'flat_plate_baseline.yaml'

# The case's own thruster definition (tcd/tdf.csv, thruster type ARG).
THRUSTER = {'d': 1.0, 've': 577.0684534784414, 'R': 208.13,
            'gamma': 1.6666666666666667, 'Te': 200, 'n': 1.0e20}

# Panel-local basis of the ISS-panel case: u = +X, v = +Y, n = +Z.
PANEL_U = np.array([1.0, 0.0, 0.0])
PANEL_V = np.array([0.0, 1.0, 0.0])
PANEL_N = np.array([0.0, 0.0, 1.0])


def baseline_mapping():
    return yaml.safe_load(BASELINE_YAML.read_text(encoding='utf-8'))


def from_mapping(data):
    return StudyConfig.from_mapping(data, source_path=str(BASELINE_YAML))


def panel_target():
    return TargetSpec.from_mapping(
        {'reference_point': [0.0, 0.0, 0.0],
         'normal': [0.0, 0.0, 1.0], 'tangent': [1.0, 0.0, 0.0]}, 'panel.stl')


def sweep_mapping(**overrides):
    """A parallel_to_normal sweep mapping with the given overrides."""
    data = {'source_axis_mode': 'parallel_to_normal',
            'source_distances': [4.0], 'n_firings': 1, 'thrusters': [1]}
    data.update(overrides)
    return data


# ---------------------------------------------------------------- model dispatch
class PlumeModelSelection(unittest.TestCase):
    """Both Cai variants are selectable; nothing else is."""

    def test_registry_holds_exactly_the_two_collisionless_models(self):
        self.assertEqual(sorted(PLUME_MODELS),
                         ['CollisionlessGasKinetics', 'SimplifiedGasKinetics'])

    def test_supported_names_resolve(self):
        for name in PLUME_MODELS:
            with self.subTest(model=name):
                self.assertEqual(resolve_model_name(name), name)

    def test_unknown_name_is_rejected_never_defaulted(self):
        for name in ('DSMC', 'Collisionless', 'simplifiedgaskinetics', ''):
            with self.subTest(model=name):
                with pytest.raises(PlumeModelError):
                    resolve_model_name(name)

    def test_omitted_name_takes_the_historical_default(self):
        self.assertEqual(resolve_model_name(None), DEFAULT_PLUME_MODEL)
        self.assertEqual(DEFAULT_PLUME_MODEL, 'SimplifiedGasKinetics')

    def test_kinetics_keys_round_trip(self):
        for key, name in (('Simplified', 'SimplifiedGasKinetics'),
                          ('Collisionless', 'CollisionlessGasKinetics')):
            with self.subTest(key=key):
                self.assertEqual(model_name_for_kinetics(key), name)
                self.assertEqual(kinetics_key_for(name), key)

    def test_disabled_and_unknown_kinetics_keys_are_reported(self):
        with pytest.raises(PlumeModelError) as excinfo:
            model_name_for_kinetics('None')
        self.assertIn('disables', str(excinfo.value))
        with pytest.raises(PlumeModelError):
            model_name_for_kinetics('DSMC')

    def test_study_configuration_accepts_both_and_rejects_others(self):
        for name in PLUME_MODELS:
            with self.subTest(model=name):
                data = baseline_mapping()
                data['plume_model']['name'] = name
                self.assertEqual(from_mapping(data).plume_model, name)

        data = baseline_mapping()
        data['plume_model']['name'] = 'BoltzmannGasKinetics'
        with pytest.raises(StudyConfigError):
            from_mapping(data)


class CommonLocalFieldState(unittest.TestCase):
    """Both models reduce to one field-state structure the GSI consumes."""

    FIELDS = ('number_density', 'mass_density', 'axial_velocity',
              'radial_velocity', 'velocity_magnitude', 'temperature',
              'speed_ratio')

    def _state(self, model_name, theta):
        model = create_model(model_name, 4.0, theta, THRUSTER, 300.0, 1.0)
        return local_field_state(model)

    def test_both_models_populate_every_field(self):
        for model_name in PLUME_MODELS:
            for theta in (0.0, 0.3):
                with self.subTest(model=model_name, theta=theta):
                    state = self._state(model_name, theta)
                    for name in self.FIELDS:
                        value = getattr(state, name)
                        self.assertTrue(np.isfinite(value),
                                        f'{name} is not finite')
                        self.assertGreaterEqual(value, 0.0)
                    self.assertEqual(state.velocity,
                                     (state.axial_velocity,
                                      state.radial_velocity))
                    self.assertIn('number_density', state.to_dict())

    def test_velocity_magnitude_is_the_axial_radial_resultant(self):
        state = self._state('CollisionlessGasKinetics', 0.3)
        self.assertAlmostEqual(
            state.velocity_magnitude,
            float(np.hypot(state.axial_velocity, state.radial_velocity)),
            places=9)

    def test_centerline_is_flagged_and_purely_axial(self):
        for model_name in PLUME_MODELS:
            with self.subTest(model=model_name):
                state = self._state(model_name, 0.0)
                self.assertTrue(state.on_centerline)
                self.assertEqual(state.radial_velocity, 0.0)

    def test_models_agree_on_the_centerline_and_differ_off_it(self):
        # The exact closed forms are shared, so the centerline must match
        # bit-for-bit; off-axis the full model integrates the exit disk and
        # must therefore differ -- that difference IS the model selection.
        centerline = [self._state(name, 0.0) for name in PLUME_MODELS]
        self.assertEqual(centerline[0].to_dict(), centerline[1].to_dict())

        simplified = self._state('SimplifiedGasKinetics', 0.3)
        full = self._state('CollisionlessGasKinetics', 0.3)
        self.assertNotAlmostEqual(simplified.number_density,
                                  full.number_density, places=6)

    def test_shared_gsi_consumes_either_state_identically(self):
        # One implementation of the Maxwellian wall formulas, fed by both
        # models: no per-model duplication of pressure/shear/heat logic.
        for model_name in PLUME_MODELS:
            with self.subTest(model=model_name):
                state = self._state(model_name, 0.2)
                pressure, shear, heat = maxwellian_surface_loads(
                    state, sigma=1.0, T_w=300.0, R=THRUSTER['R'],
                    gamma=THRUSTER['gamma'], incidence=0.2)
                self.assertGreater(pressure, 0.0)
                self.assertGreater(heat, 0.0)
                self.assertTrue(np.isfinite(shear))


# --------------------------------------------------------------- pose generation
class PanelLocalPoseGeneration(unittest.TestCase):
    """parallel_to_normal translates the source without tilting the axis."""

    def test_centered_normal_incidence_pose(self):
        position, dcm = translated_pose_for(4.0, 0.0, 0.0, [0, 0, 0],
                                            PANEL_N, PANEL_U)
        np.testing.assert_allclose(position, [0.0, 0.0, 4.0], atol=1e-12)
        # First column is the thruster axis, anti-parallel to the normal.
        np.testing.assert_allclose(dcm[:, 0], -PANEL_N, atol=1e-12)
        np.testing.assert_allclose(dcm[:, 1], PANEL_V, atol=1e-12)
        np.testing.assert_allclose(dcm[:, 2], PANEL_U, atol=1e-12)

    def test_dcm_is_a_proper_orthonormal_rotation(self):
        for offset_u, offset_v in ((0.0, 0.0), (-9.0, 0.0), (5.5, 2.0)):
            with self.subTest(u=offset_u, v=offset_v):
                _, dcm = translated_pose_for(4.0, offset_u, offset_v,
                                             [0, 0, 0], PANEL_N, PANEL_U)
                np.testing.assert_allclose(dcm.T @ dcm, np.eye(3), atol=1e-12)
                self.assertAlmostEqual(float(np.linalg.det(dcm)), 1.0,
                                       places=12)

    def test_u_offset_translates_parallel_and_leaves_the_axis_alone(self):
        _, reference_dcm = translated_pose_for(4.0, 0.0, 0.0, [0, 0, 0],
                                               PANEL_N, PANEL_U)
        for offset_u in (-9.0, -4.5, 4.5, 9.0):
            with self.subTest(u=offset_u):
                position, dcm = translated_pose_for(
                    4.0, offset_u, 0.0, [0, 0, 0], PANEL_N, PANEL_U)
                # Translated along u only; stand-off unchanged.
                np.testing.assert_allclose(position,
                                           [offset_u, 0.0, 4.0], atol=1e-12)
                # The axis is still exactly -n: NOT re-aimed at the center.
                np.testing.assert_allclose(dcm[:, 0], -PANEL_N, atol=1e-12)
                np.testing.assert_allclose(dcm, reference_dcm, atol=1e-12)

    def test_v_offset_translates_along_the_transverse_axis(self):
        position, dcm = translated_pose_for(6.0, 0.0, 3.0, [0, 0, 0],
                                            PANEL_N, PANEL_U)
        np.testing.assert_allclose(position, [0.0, 3.0, 6.0], atol=1e-12)
        np.testing.assert_allclose(dcm[:, 0], -PANEL_N, atol=1e-12)

    def test_zero_offset_reproduces_the_aim_at_reference_head_on_pose(self):
        # The new mode EXTENDS the old convention rather than redefining it.
        aimed_position, aimed_dcm = pose_for(0.0, 4.0, [0, 0, 0],
                                             PANEL_N, PANEL_U)
        position, dcm = translated_pose_for(4.0, 0.0, 0.0, [0, 0, 0],
                                            PANEL_N, PANEL_U)
        np.testing.assert_allclose(position, aimed_position, atol=1e-12)
        np.testing.assert_allclose(dcm, aimed_dcm, atol=1e-12)

    def test_aimed_mode_is_untouched_by_the_new_dispatcher(self):
        pose = SweepPose(plate_angle_deg=30.0, source_distance=5.0,
                         axis_mode='aim_at_reference')
        target = panel_target()
        position, dcm = pose_for_sweep_pose(pose, target)
        expected = pose_for(30.0, 5.0, target.reference_point, target.normal,
                            target.tangent)
        np.testing.assert_allclose(position, expected[0], atol=1e-12)
        np.testing.assert_allclose(dcm, expected[1], atol=1e-12)

    def test_generated_firings_carry_the_offsets(self):
        sweep = SweepSpec.from_mapping(
            sweep_mapping(source_offsets_u=[7.0], source_offsets_v=[-2.0]))
        pose = sweep.sweep_poses[0]
        firings = build_case_firings(sweep, panel_target(),
                                     pose.plate_angle_deg,
                                     pose.source_distance, pose=pose)
        self.assertEqual(len(firings), 1)
        self.assertEqual(firings[0].source_offset_u, 7.0)
        self.assertEqual(firings[0].source_offset_v, -2.0)
        self.assertEqual(firings[0].source_axis_mode, 'parallel_to_normal')
        np.testing.assert_allclose(firings[0].position, [7.0, -2.0, 4.0],
                                   atol=1e-12)


class OffsetSweepEnumeration(unittest.TestCase):
    """Case count and order for the offset grid."""

    def test_offsets_default_to_zero(self):
        sweep = SweepSpec.from_mapping({'source_distances': [4.0]})
        self.assertEqual(sweep.source_offsets_u, (0.0,))
        self.assertEqual(sweep.source_offsets_v, (0.0,))
        self.assertEqual(sweep.source_axis_mode, 'aim_at_reference')
        for pose in sweep.sweep_poses:
            self.assertEqual(pose.source_offset_u, 0.0)
            self.assertEqual(pose.source_offset_v, 0.0)

    def test_case_count_is_the_product_of_the_three_axes(self):
        distances = [2.0, 4.0, 6.0]
        offsets_u = [-9.0, -4.5, 0.0, 4.5, 9.0]
        offsets_v = [0.0, 3.0]
        sweep = SweepSpec.from_mapping(sweep_mapping(
            source_distances=distances, source_offsets_u=offsets_u,
            source_offsets_v=offsets_v))
        self.assertEqual(len(sweep.sweep_poses),
                         len(distances) * len(offsets_u) * len(offsets_v))
        self.assertEqual(sweep.total_firings, len(sweep.sweep_poses))

    def test_pose_order_is_distance_then_u_then_v(self):
        sweep = SweepSpec.from_mapping(sweep_mapping(
            source_distances=[2.0, 4.0], source_offsets_u=[-1.0, 1.0],
            source_offsets_v=[0.0, 5.0]))
        self.assertEqual(
            [(p.source_distance, p.source_offset_u, p.source_offset_v)
             for p in sweep.sweep_poses],
            [(2.0, -1.0, 0.0), (2.0, -1.0, 5.0),
             (2.0, 1.0, 0.0), (2.0, 1.0, 5.0),
             (4.0, -1.0, 0.0), (4.0, -1.0, 5.0),
             (4.0, 1.0, 0.0), (4.0, 1.0, 5.0)])

    def test_default_offsets_preserve_the_historical_pose_order(self):
        sweep = SweepSpec.from_mapping(
            {'plate_angles_deg': [-10.0, 10.0],
             'source_distances': [2.0, 4.0]})
        self.assertEqual(sweep.poses,
                         ((-10.0, 2.0), (10.0, 2.0),
                          (-10.0, 4.0), (10.0, 4.0)))

    def test_non_finite_offsets_are_rejected(self):
        for bad in ([float('nan')], [float('inf')], [0.0, float('-inf')]):
            with self.subTest(offsets=bad):
                with pytest.raises(StudyConfigError):
                    SweepSpec.from_mapping(sweep_mapping(
                        source_offsets_u=bad))

    def test_empty_offset_list_is_rejected(self):
        with pytest.raises(StudyConfigError):
            SweepSpec.from_mapping(sweep_mapping(source_offsets_u=[]))

    def test_unknown_axis_mode_is_rejected(self):
        with pytest.raises(StudyConfigError) as excinfo:
            SweepSpec.from_mapping({'source_distances': [4.0],
                                    'source_axis_mode': 'follow_the_plume'})
        self.assertIn('source_axis_mode', str(excinfo.value))

    def test_incompatible_pose_definitions_are_never_silently_combined(self):
        # Offsets while aiming at the reference point: the offset has no
        # meaning when the axis is re-aimed every time.
        with pytest.raises(StudyConfigError) as excinfo:
            SweepSpec.from_mapping({'source_distances': [4.0],
                                    'source_offsets_u': [3.0]})
        self.assertIn('parallel_to_normal', str(excinfo.value))

        # An approach angle with a fixed axis: the angle cannot be realized.
        with pytest.raises(StudyConfigError) as excinfo:
            SweepSpec.from_mapping(sweep_mapping(
                plate_angles_deg=[0.0, 30.0]))
        self.assertIn('plate_angles_deg', str(excinfo.value))


# ------------------------------------------------------------ Knudsen metadata
class DerivedKnudsenMetadata(unittest.TestCase):
    """Kn is computed correctly, and is metadata only."""

    def test_source_distance_reference_length(self):
        spec = KnudsenSpec.from_mapping({
            'mean_free_path_m': 1.0, 'reference_length': 'source_distance'})
        self.assertEqual(spec.reference_mode, 'source_distance')
        self.assertEqual(spec.definition, 'lambda_over_source_distance')
        for distance, expected in ((2.0, 0.5), (4.0, 0.25), (10.0, 0.1)):
            with self.subTest(L=distance):
                self.assertAlmostEqual(spec.knudsen_number(distance),
                                       expected, places=12)
                self.assertEqual(spec.reference_length_for(distance), distance)

    def test_explicit_reference_length(self):
        spec = KnudsenSpec.from_mapping({
            'mean_free_path_m': 0.1, 'reference_length_m': 1.0,
            'definition': 'lambda_over_nozzle_diameter'})
        self.assertEqual(spec.reference_mode, 'explicit')
        self.assertEqual(spec.definition, 'lambda_over_nozzle_diameter')
        # Fixed length: the same Kn whatever the case's stand-off is.
        for distance in (2.0, 4.0, 100.0):
            with self.subTest(L=distance):
                self.assertAlmostEqual(spec.knudsen_number(distance), 0.1,
                                       places=12)
                self.assertEqual(spec.reference_length_for(distance), 1.0)

    def test_the_documented_kn_labels_are_reproducible(self):
        # One study per label, nozzle diameter D = 1 m as the reference.
        for mean_free_path, label in ((100.0, 100.0), (10.0, 10.0),
                                      (1.0, 1.0), (0.1, 0.1), (0.01, 0.01)):
            with self.subTest(Kn=label):
                spec = KnudsenSpec.from_mapping({
                    'mean_free_path_m': mean_free_path,
                    'reference_length_m': 1.0})
                self.assertAlmostEqual(spec.knudsen_number(4.0), label,
                                       places=12)

    def test_omitted_block_yields_no_spec_and_no_fields(self):
        self.assertIsNone(KnudsenSpec.from_mapping(None))
        self.assertIsNone(KnudsenSpec.from_mapping({}))
        config = from_mapping(baseline_mapping())
        self.assertIsNone(config.knudsen)
        self.assertNotIn('knudsen', config.provenance())

    def test_missing_mean_free_path_is_refused_never_inferred(self):
        with pytest.raises(StudyConfigError) as excinfo:
            KnudsenSpec.from_mapping({'reference_length': 'source_distance'})
        message = str(excinfo.value)
        self.assertIn('mean_free_path_m', message)
        self.assertIn('never infers', message)

    def test_non_positive_or_non_finite_mean_free_path_is_rejected(self):
        for value in (0.0, -1.0, float('nan'), float('inf'), 'thin'):
            with self.subTest(mean_free_path=value):
                with pytest.raises(StudyConfigError):
                    KnudsenSpec.from_mapping({
                        'mean_free_path_m': value,
                        'reference_length': 'source_distance'})

    def test_exactly_one_reference_length_mode_is_required(self):
        # Neither.
        with pytest.raises(StudyConfigError) as excinfo:
            KnudsenSpec.from_mapping({'mean_free_path_m': 1.0})
        self.assertIn('EXACTLY ONE', str(excinfo.value))
        # Both.
        with pytest.raises(StudyConfigError) as excinfo:
            KnudsenSpec.from_mapping({'mean_free_path_m': 1.0,
                                      'reference_length': 'source_distance',
                                      'reference_length_m': 1.0})
        self.assertIn('EXACTLY ONE', str(excinfo.value))

    def test_unknown_symbolic_reference_length_is_rejected(self):
        with pytest.raises(StudyConfigError) as excinfo:
            KnudsenSpec.from_mapping({'mean_free_path_m': 1.0,
                                      'reference_length': 'nozzle_diameter'})
        self.assertIn('reference_length_m', str(excinfo.value))

    def test_invalid_explicit_reference_length_is_rejected(self):
        for value in (0.0, -2.0, float('inf'), 'wide'):
            with self.subTest(reference_length_m=value):
                with pytest.raises(StudyConfigError):
                    KnudsenSpec.from_mapping({'mean_free_path_m': 1.0,
                                              'reference_length_m': value})

    def test_provenance_records_kn_as_derived_metadata(self):
        data = baseline_mapping()
        data['knudsen'] = {'mean_free_path_m': 1.0,
                           'reference_length': 'source_distance'}
        provenance = from_mapping(data).provenance()
        self.assertIn('knudsen', provenance)
        self.assertIn('derived metadata only',
                      provenance['knudsen']['role'])


# ---------------------------------------------------- panel-local projections
class PanelLocalCoordinates(unittest.TestCase):

    def test_known_coordinates_are_recovered(self):
        points = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0],
                           [0.0, -2.0, 0.0], [-1.5, 4.0, 0.0]])
        local_u, local_v = panel_local_coordinates(points, [0, 0, 0],
                                                   PANEL_U, PANEL_V)
        np.testing.assert_allclose(local_u, [0.0, 3.0, 0.0, -1.5], atol=1e-12)
        np.testing.assert_allclose(local_v, [0.0, 0.0, -2.0, 4.0], atol=1e-12)

    def test_coordinates_are_measured_from_the_reference_point(self):
        points = np.array([[5.0, 7.0, 0.0]])
        local_u, local_v = panel_local_coordinates(points, [2.0, 3.0, 0.0],
                                                   PANEL_U, PANEL_V)
        np.testing.assert_allclose(local_u, [3.0], atol=1e-12)
        np.testing.assert_allclose(local_v, [4.0], atol=1e-12)

    def test_a_rotated_basis_still_recovers_the_panel_coordinates(self):
        # A panel in the X-Z plane: u = +X, n = -Y, v = n x u = +Z.
        u_hat = np.array([1.0, 0.0, 0.0])
        n_hat = np.array([0.0, -1.0, 0.0])
        v_hat = np.cross(n_hat, u_hat)
        np.testing.assert_allclose(v_hat, [0.0, 0.0, 1.0], atol=1e-12)
        local_u, local_v = panel_local_coordinates(
            np.array([[2.0, 0.0, -3.0]]), [0, 0, 0], u_hat, v_hat)
        np.testing.assert_allclose(local_u, [2.0], atol=1e-12)
        np.testing.assert_allclose(local_v, [-3.0], atol=1e-12)


def _loads(centroids, pressures, *, areas=None, reference=(0.0, 0.0, 0.0),
           source=(0.0, 0.0, 4.0)):
    """Integrate a synthetic pressure-only load on a +Z-facing panel."""
    centroids = np.asarray(centroids, dtype=float)
    n = len(centroids)
    areas = np.ones(n) if areas is None else np.asarray(areas, dtype=float)
    pressures = np.asarray(pressures, dtype=float)
    return integrate_component_loads(
        component_name='panel', face_indices=np.arange(n),
        centroids=centroids, unit_normals=np.tile(PANEL_N, (n, 1)),
        areas=areas, pressures=pressures, shear_stresses=np.zeros(n),
        heat_fluxes=np.zeros(n), strikes=(pressures != 0.0).astype(float),
        moment_reference_point=reference, source_position=source)


class PanelLoadProjection(unittest.TestCase):
    """Hand-computed normal force, moments and center of pressure."""

    def test_normal_force_is_positive_into_the_panel(self):
        # 10 Pa over 1 m^2 on a +Z-facing face: the plume pushes along -Z,
        # so the global force is -10 Z and the normal force is +10 N.
        loads = _loads([[0.0, 0.0, 0.0]], [10.0])
        np.testing.assert_allclose(loads.force, [0.0, 0.0, -10.0], atol=1e-12)
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(panel.normal_force, 10.0, places=12)
        self.assertAlmostEqual(panel.local_force_u, 0.0, places=12)
        self.assertAlmostEqual(panel.local_force_v, 0.0, places=12)

    def test_symmetric_centered_load_gives_no_moment_about_the_center(self):
        loads = _loads([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                        [0.0, -2.0, 0.0], [0.0, 2.0, 0.0]],
                       [7.0, 7.0, 7.0, 7.0])
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(panel.local_moment_u, 0.0, places=12)
        self.assertAlmostEqual(panel.local_moment_v, 0.0, places=12)
        self.assertAlmostEqual(panel.local_moment_n, 0.0, places=12)
        self.assertAlmostEqual(panel.center_of_pressure_u, 0.0, places=12)
        self.assertAlmostEqual(panel.center_of_pressure_v, 0.0, places=12)
        self.assertAlmostEqual(panel.normal_force, 28.0, places=12)

    def test_off_center_load_gives_the_right_moment_sign_and_magnitude(self):
        # One 10 Pa face of 1 m^2 at u = +2: F = -10 Z, arm = +2 X, so
        # M = (2 X) x (-10 Z) = +20 Y = +20 v.  A source displaced toward
        # +u therefore gives a POSITIVE local_moment_v.
        loads = _loads([[2.0, 0.0, 0.0]], [10.0])
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(panel.local_moment_v, 20.0, places=12)
        self.assertAlmostEqual(panel.local_moment_u, 0.0, places=12)
        self.assertAlmostEqual(panel.center_of_pressure_u, 2.0, places=12)

        # Mirrored offset: equal magnitude, opposite sign.
        mirrored = project_to_panel_frame(_loads([[-2.0, 0.0, 0.0]], [10.0]),
                                          PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(mirrored.local_moment_v, -20.0, places=12)
        self.assertAlmostEqual(mirrored.center_of_pressure_u, -2.0, places=12)

    def test_transverse_offset_loads_the_u_moment(self):
        # A face at v = +3: M = (3 Y) x (-10 Z) = -30 X = -30 u.
        loads = _loads([[0.0, 3.0, 0.0]], [10.0])
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(panel.local_moment_u, -30.0, places=12)
        self.assertAlmostEqual(panel.local_moment_v, 0.0, places=12)
        self.assertAlmostEqual(panel.center_of_pressure_v, 3.0, places=12)

    def test_center_of_pressure_is_measured_from_the_moment_reference(self):
        loads = _loads([[5.0, 0.0, 0.0]], [10.0], reference=(3.0, 0.0, 0.0))
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertAlmostEqual(panel.center_of_pressure_u, 2.0, places=12)

    def test_unavailable_center_of_pressure_projects_to_none(self):
        loads = _loads([[1.0, 0.0, 0.0]], [0.0])
        panel = project_to_panel_frame(loads, PANEL_U, PANEL_V, PANEL_N)
        self.assertIsNone(panel.center_of_pressure_u)
        self.assertIsNone(panel.center_of_pressure_v)
        self.assertAlmostEqual(panel.normal_force, 0.0, places=12)


# ------------------------------------------------------ distribution export
class SurfaceDistributionExport(unittest.TestCase):

    def _mesh(self, n=6):
        centroids = np.column_stack([np.linspace(-5.0, 5.0, n),
                                     np.linspace(-2.0, 2.0, n),
                                     np.zeros(n)])
        return centroids, np.full(n, 0.5)

    def test_rows_cover_every_face_with_the_required_columns(self):
        centroids, areas = self._mesh()
        n = len(centroids)
        rows = distribution_rows(
            np.arange(n), centroids, areas, np.linspace(1.0, 2.0, n),
            np.linspace(0.1, 0.2, n), np.linspace(10.0, 20.0, n),
            np.ones(n), reference_point=[0, 0, 0], u_hat=PANEL_U,
            v_hat=PANEL_V)
        self.assertEqual(len(rows), n)
        for row in rows:
            self.assertEqual(set(row), set(DISTRIBUTION_COLUMNS))

    def test_local_coordinates_match_the_panel_basis(self):
        centroids, areas = self._mesh()
        n = len(centroids)
        rows = distribution_rows(
            np.arange(n), centroids, areas, np.zeros(n), np.zeros(n),
            np.zeros(n), None, reference_point=[0, 0, 0], u_hat=PANEL_U,
            v_hat=PANEL_V)
        np.testing.assert_allclose([row['local_u'] for row in rows],
                                   centroids[:, 0], atol=1e-12)
        np.testing.assert_allclose([row['local_v'] for row in rows],
                                   centroids[:, 1], atol=1e-12)

    def test_component_subset_preserves_the_original_face_indices(self):
        centroids, areas = self._mesh()
        n = len(centroids)
        rows = distribution_rows(
            [1, 4], centroids, areas, np.arange(n, dtype=float),
            np.zeros(n), np.zeros(n), None, reference_point=[0, 0, 0],
            u_hat=PANEL_U, v_hat=PANEL_V)
        self.assertEqual([row['face_index'] for row in rows], [1, 4])
        self.assertEqual([row['pressure'] for row in rows], [1.0, 4.0])

    def test_values_are_written_without_interpolation(self):
        centroids, areas = self._mesh()
        n = len(centroids)
        pressures = np.linspace(3.0, 9.0, n)
        rows = distribution_rows(
            np.arange(n), centroids, areas, pressures, np.zeros(n),
            np.zeros(n), None, reference_point=[0, 0, 0], u_hat=PANEL_U,
            v_hat=PANEL_V)
        with tempfile.TemporaryDirectory() as tmp:
            path = write_surface_distribution(
                os.path.join(tmp, 'dist.csv'), rows,
                {'study_name': 'unit', 'plume_model': 'SimplifiedGasKinetics'})
            with open(path, encoding='utf-8', newline='') as handle:
                written = list(csv.DictReader(handle))
            self.assertEqual(len(written), n)
            self.assertEqual(list(written[0]), list(DISTRIBUTION_COLUMNS))
            np.testing.assert_allclose(
                [float(row['pressure']) for row in written], pressures,
                rtol=0, atol=0)

            sidecar = os.path.splitext(path)[0] + '.meta.json'
            self.assertTrue(os.path.isfile(sidecar))
            import json
            document = json.loads(Path(sidecar).read_text(encoding='utf-8'))
            self.assertEqual(document['units']['pressure'], 'Pa')
            self.assertEqual(document['units']['local_u'], 'm')
            self.assertEqual(document['n_faces'], n)
            self.assertEqual(document['plume_model'], 'SimplifiedGasKinetics')

    def test_empty_distribution_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with pytest.raises(ValueError):
                write_surface_distribution(os.path.join(tmp, 'x.csv'), [])


# --------------------------------------------------------- result serialization
def _case_result(**overrides):
    fields = dict(
        study_name='s', case_id='case000', component='panel', firing_id=1,
        geometry_id='panel.stl', mesh_faces=2, component_faces=2,
        component_area=1.0, coordinate_system='global', units={'length': 'm'},
        plume_source_position=[0.0, 0.0, 4.0],
        plume_source_orientation=[1.0] + [0.0] * 8,
        target_normal=[0.0, 0.0, 1.0], target_tangent=[1.0, 0.0, 0.0],
        target_reference_point=[0.0, 0.0, 0.0], plate_angle_deg=0.0,
        source_distance=4.0, firing_duration_s=1.0, thrusters=[1],
        plume_model='CollisionlessGasKinetics', plume_model_parameters={},
        pressure_force=[0.0, 0.0, -1.0], shear_force=[0.0, 0.0, 0.0],
        force=[0.0, 0.0, -1.0], force_magnitude=1.0,
        moment_reference_point=[0.0, 0.0, 0.0],
        pressure_moment=[0.0, 0.0, 0.0], shear_moment=[0.0, 0.0, 0.0],
        moment=[0.0, 2.0, 0.0], moment_magnitude=2.0,
        center_of_pressure=[2.0, 0.0, 0.0], center_of_pressure_status='ok',
        residual_couple=0.0, pressure_weighted_centroid=[2.0, 0.0, 0.0],
        max_pressure=0.3, max_shear_stress=0.03, max_heat_flux=50.0,
        total_heat_load=100.0, affected_area=0.5, struck_faces=1)
    fields.update(overrides)
    return CaseResult(**fields)


class ResultSchemaExtensions(unittest.TestCase):

    def test_new_fields_default_safely(self):
        # Every addition is optional, so an older-style construction works.
        case = _case_result()
        self.assertEqual(case.source_offset_u, 0.0)
        self.assertEqual(case.source_offset_v, 0.0)
        self.assertEqual(case.source_axis_mode, 'aim_at_reference')
        self.assertIsNone(case.normal_force)
        self.assertIsNone(case.knudsen_number)
        self.assertIsNone(case.surface_distribution_path)

    def test_model_variant_is_derived_from_the_model_name(self):
        self.assertEqual(_case_result().model_variant, 'Collisionless')
        self.assertEqual(
            _case_result(plume_model='SimplifiedGasKinetics').model_variant,
            'Simplified')

    def test_csv_row_carries_the_new_metadata_flat(self):
        case = _case_result(
            source_offset_u=5.5, source_offset_v=-1.0,
            source_axis_mode='parallel_to_normal', normal_force=2.5,
            local_moment_u=0.0, local_moment_v=13.75, local_moment_n=0.0,
            center_of_pressure_u=5.4, center_of_pressure_v=-1.0,
            knudsen_number=0.25, mean_free_path=1.0,
            knudsen_reference_length=4.0,
            knudsen_definition='lambda_over_source_distance',
            surface_distribution_path='/tmp/dist.csv')
        row = case.to_row()
        self.assertEqual(row['source_offset_u'], 5.5)
        self.assertEqual(row['source_offset_v'], -1.0)
        self.assertEqual(row['source_axis_mode'], 'parallel_to_normal')
        self.assertEqual(row['model_variant'], 'Collisionless')
        self.assertEqual(row['normal_force'], 2.5)
        self.assertEqual(row['local_moment_v'], 13.75)
        self.assertEqual(row['center_of_pressure_u'], 5.4)
        self.assertEqual(row['knudsen_number'], 0.25)
        self.assertEqual(row['knudsen_definition'],
                         'lambda_over_source_distance')
        self.assertEqual(row['surface_distribution_path'], '/tmp/dist.csv')
        # Flat: no nested containers or arrays in a CSV row.
        for key, value in row.items():
            self.assertNotIsInstance(value, (list, dict, tuple, np.ndarray),
                                     f'column {key} is not flat')

    def test_absent_optional_values_leave_empty_columns(self):
        row = _case_result().to_row()
        for column in ('normal_force', 'local_moment_v', 'knudsen_number',
                       'mean_free_path', 'knudsen_reference_length',
                       'knudsen_definition', 'surface_distribution_path'):
            self.assertEqual(row[column], '', f'{column} should be empty')

    def test_json_form_keeps_model_and_knudsen_provenance(self):
        case = _case_result(knudsen_number=0.25, mean_free_path=1.0,
                            knudsen_definition='lambda_over_source_distance')
        document = case.to_dict()
        self.assertEqual(document['plume_model'], 'CollisionlessGasKinetics')
        self.assertEqual(document['knudsen_number'], 0.25)
        self.assertEqual(document['knudsen_definition'],
                         'lambda_over_source_distance')

    def test_quantity_exposes_the_new_comparable_scalars(self):
        case = _case_result(normal_force=2.5, local_moment_v=13.75,
                            center_of_pressure_u=5.4, knudsen_number=0.25)
        self.assertEqual(case.quantity('normal_force'), 2.5)
        self.assertEqual(case.quantity('local_moment_v'), 13.75)
        self.assertEqual(case.quantity('center_of_pressure_u'), 5.4)
        self.assertEqual(case.quantity('knudsen_number'), 0.25)
        # Unavailable quantities stay None rather than becoming zero.
        self.assertIsNone(_case_result().quantity('normal_force'))
        self.assertIsNone(case.quantity('not_a_quantity'))


# ------------------------------------------------- committed ISS configurations
class CommittedISSPanelConfigurations(unittest.TestCase):

    def _config(self, name):
        return StudyConfig.from_yaml(ISS_CASE / 'study' / name)

    def test_baseline_configurations_differ_only_in_the_model(self):
        simplified = self._config('iss_panel_baseline_simplified.yaml')
        full = self._config('iss_panel_baseline_full_cai.yaml')
        self.assertEqual(simplified.plume_model, 'SimplifiedGasKinetics')
        self.assertEqual(full.plume_model, 'CollisionlessGasKinetics')
        # Same geometry, same pose, same loads definition.
        self.assertEqual(simplified.sweep.sweep_poses,
                         full.sweep.sweep_poses)
        np.testing.assert_allclose(simplified.target.reference_point,
                                   full.target.reference_point)
        self.assertEqual(simplified.loads.normalization.reference_area,
                         full.loads.normalization.reference_area)

    def test_baseline_is_a_centered_parallel_to_normal_pose(self):
        config = self._config('iss_panel_baseline_simplified.yaml')
        self.assertEqual(config.sweep.source_axis_mode, 'parallel_to_normal')
        self.assertEqual(config.n_cases, 1)
        pose = config.sweep.sweep_poses[0]
        self.assertEqual((pose.source_distance, pose.source_offset_u,
                          pose.source_offset_v), (4.0, 0.0, 0.0))
        position, dcm = pose_for_sweep_pose(pose, config.target)
        np.testing.assert_allclose(position, [0.0, 0.0, 4.0], atol=1e-12)
        np.testing.assert_allclose(dcm[:, 0], [0.0, 0.0, -1.0], atol=1e-12)

    def test_panel_basis_is_the_documented_one(self):
        config = self._config('iss_panel_baseline_simplified.yaml')
        u_hat, v_hat, n_hat = config.target.local_basis()
        np.testing.assert_allclose(u_hat, PANEL_U, atol=1e-12)
        np.testing.assert_allclose(v_hat, PANEL_V, atol=1e-12)
        np.testing.assert_allclose(n_hat, PANEL_N, atol=1e-12)
        np.testing.assert_allclose(np.cross(u_hat, v_hat), n_hat, atol=1e-12)

    def test_sweep_enumerates_distances_times_offsets(self):
        config = self._config('iss_panel_offset_distance_sweep.yaml')
        self.assertEqual(config.sweep.mode, 'per_case')
        self.assertEqual(len(config.sweep.source_distances), 3)
        self.assertEqual(len(config.sweep.source_offsets_u), 5)
        self.assertEqual(len(config.sweep.source_offsets_v), 1)
        self.assertEqual(config.n_cases, 15)
        self.assertEqual(len(config.sweep.sweep_poses), 15)

    def test_sweep_configures_distributions_plots_and_knudsen(self):
        config = self._config('iss_panel_offset_distance_sweep.yaml')
        self.assertTrue(config.output.write_surface_distribution)
        self.assertTrue(config.output.write_plots)
        self.assertTrue(config.output.write_distribution_plots)
        self.assertIsNotNone(config.knudsen)
        self.assertEqual(config.knudsen.reference_mode, 'source_distance')
        self.assertAlmostEqual(config.knudsen.knudsen_number(2.0), 0.5)
        self.assertAlmostEqual(config.knudsen.knudsen_number(4.0), 0.25)

    def test_offsets_stay_within_the_panel(self):
        config = self._config('iss_panel_offset_distance_sweep.yaml')
        for offset in config.sweep.source_offsets_u:
            self.assertLessEqual(abs(offset), 11.0)
        for offset in config.sweep.source_offsets_v:
            self.assertLessEqual(abs(offset), 6.0)


if __name__ == '__main__':
    unittest.main()
