# ========================
# PyRPOD: tests/mdao/mdao_integration_test_05.py
# ========================
# The committed ISS-representative solar-panel case run end to end through
# the package-level trade-study API:
#
#     study = TradeStudy.from_config('.../iss_panel_baseline_simplified.yaml')
#     results = study.run()
#
# The case is an idealized 22 m x 12 m flat panel (case/plume/iss_panel_thesis)
# struck by the Cai 2016 argon round jet, with the source TRANSLATED parallel
# to the panel (source_axis_mode: parallel_to_normal). Checked here:
#
#   * the committed baseline runs and produces the panel-local resultants,
#     the derived Knudsen metadata and the distribution export;
#   * selecting CollisionlessGasKinetics changes the ANSWER, not just the
#     recorded metadata -- the one property that proves model dispatch
#     reaches the calculation;
#   * the integrated pressure force recovered from the exported per-face
#     values matches the CaseResult pressure force;
#   * a centered source produces a symmetric load: no moment about the panel
#     center and a center of pressure at the origin; an offset source moves
#     both, with the documented sign;
#   * the offset sweep enumerates n_distances x n_u_offsets x n_v_offsets
#     cases with stable identifiers, and its CSV/JSON carry the new columns;
#   * plot generation completes headless.
#
# Small meshes and small sweeps: the offset-sweep tests build their own
# coarse panel and a 2 x 3 grid rather than running the committed 15-case
# sweep on 2112 faces. The BASELINE runs on the committed mesh unchanged.
#
# No DSMC is involved anywhere: both models under test are collisionless and
# Kn is asserted to be metadata only.
#
# Run:  python -m pytest mdao/mdao_integration_test_05.py   (from tests/)

import copy
import csv
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

from pyrpod.mdao.TradeStudy import TradeStudy
from pyrpod.mdao.study_config import StudyConfig

_TESTS_DIR = Path(__file__).resolve().parents[1]
CASE_DIR = _TESTS_DIR.parent / 'case' / 'plume' / 'iss_panel_thesis'
STUDY_DIR = CASE_DIR / 'study'
BASELINE_YAML = STUDY_DIR / 'iss_panel_baseline_simplified.yaml'
FULL_CAI_YAML = STUDY_DIR / 'iss_panel_baseline_full_cai.yaml'
SWEEP_YAML = STUDY_DIR / 'iss_panel_offset_distance_sweep.yaml'

PANEL_LENGTH_U = 22.0
PANEL_WIDTH_V = 12.0
PANEL_AREA = PANEL_LENGTH_U * PANEL_WIDTH_V

PANEL_U = np.array([1.0, 0.0, 0.0])
PANEL_V = np.array([0.0, 1.0, 0.0])
PANEL_N = np.array([0.0, 0.0, 1.0])


def run_config(mapping, output_dir, name):
    """Write a modified configuration into the case's study dir and run it.

    The file must live beside the committed ones so its relative
    ``case_dir: ..`` still resolves to the case; it is removed afterwards.
    """
    path = STUDY_DIR / f'_tmp_{name}.yaml'
    path.write_text(yaml.safe_dump(mapping, sort_keys=False),
                    encoding='utf-8')
    try:
        study = TradeStudy.from_config(path, output_dir=output_dir)
        return study, study.run()
    finally:
        path.unlink(missing_ok=True)


def coarse_panel_mapping(**sweep_overrides):
    """The baseline configuration on a coarse mesh, for the sweep tests."""
    mapping = yaml.safe_load(BASELINE_YAML.read_text(encoding='utf-8'))
    mapping['study']['name'] = 'iss_panel_coarse'
    mapping['target']['geometry_id'] = 'iss_panel_coarse.stl'
    mapping['sweep'].update(sweep_overrides)
    return mapping


class ISSPanelBaseline(unittest.TestCase):
    """The committed simplified baseline, on the committed mesh."""

    @classmethod
    def setUpClass(cls):
        cls.output_dir = tempfile.mkdtemp(prefix='pyrpod_iss_baseline_')
        cls.study = TradeStudy.from_config(BASELINE_YAML,
                                           output_dir=cls.output_dir)
        cls.results = cls.study.run()
        cls.case = cls.results.cases[0]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.output_dir, ignore_errors=True)

    # ------------------------------------------------------------ structure
    def test_one_case_one_component_one_firing(self):
        self.assertEqual(len(self.results), 1)
        self.assertEqual(self.case.study_name, 'iss_panel_baseline_simplified')
        self.assertEqual(self.case.component, 'panel')
        self.assertEqual(self.case.firing_id, 1)
        self.assertEqual(self.case.plume_model, 'SimplifiedGasKinetics')

    def test_target_mesh_is_the_committed_panel(self):
        self.assertAlmostEqual(self.case.component_area, PANEL_AREA, places=6)
        self.assertEqual(self.case.mesh_faces, self.case.component_faces)

    def test_pose_is_centered_and_parallel_to_the_panel_normal(self):
        self.assertEqual(self.case.source_axis_mode, 'parallel_to_normal')
        self.assertEqual(self.case.source_offset_u, 0.0)
        self.assertEqual(self.case.source_offset_v, 0.0)
        np.testing.assert_allclose(self.case.plume_source_position,
                                   [0.0, 0.0, 4.0], atol=1e-9)
        # First DCM column is the thruster axis: exactly -n, not re-aimed.
        dcm = np.asarray(self.case.plume_source_orientation).reshape(3, 3)
        np.testing.assert_allclose(dcm[:, 0], -PANEL_N, atol=1e-9)

    def test_case_id_names_the_model_distance_and_offsets(self):
        self.assertEqual(self.case.case_id,
                         'case000_modelSimplified_L4_u0_v0')

    # -------------------------------------------------------- panel-local loads
    def test_normal_force_is_positive_into_the_panel(self):
        # The plume pushes along -n, so the global force is -Z and the
        # reported normal force is positive.
        self.assertLess(self.case.force[2], 0.0)
        self.assertGreater(self.case.normal_force, 0.0)
        self.assertAlmostEqual(self.case.normal_force,
                               -float(self.case.force[2]), places=9)

    def test_centered_source_gives_no_moment_about_the_panel_center(self):
        # The panel is symmetric about its center and the source is on the
        # normal through it, so every local moment must vanish.
        for name in ('local_moment_u', 'local_moment_v', 'local_moment_n'):
            with self.subTest(component=name):
                self.assertAlmostEqual(getattr(self.case, name), 0.0,
                                       delta=1e-9 * self.case.normal_force
                                       * PANEL_LENGTH_U)

    def test_centered_source_puts_the_center_of_pressure_at_the_origin(self):
        self.assertEqual(self.case.center_of_pressure_status, 'ok')
        self.assertAlmostEqual(self.case.center_of_pressure_u, 0.0, places=6)
        self.assertAlmostEqual(self.case.center_of_pressure_v, 0.0, places=6)

    def test_surface_field_peaks_are_positive_and_finite(self):
        self.assertGreater(self.case.max_pressure, 0.0)
        self.assertGreater(self.case.max_heat_flux, 0.0)
        self.assertGreater(self.case.affected_area, 0.0)
        self.assertLessEqual(self.case.affected_area, PANEL_AREA + 1e-9)

    # ------------------------------------------------------ Knudsen metadata
    def test_knudsen_is_derived_from_the_configured_mean_free_path(self):
        self.assertAlmostEqual(self.case.mean_free_path, 1.0, places=12)
        self.assertAlmostEqual(self.case.knudsen_reference_length, 4.0,
                               places=12)
        self.assertAlmostEqual(self.case.knudsen_number, 0.25, places=12)
        self.assertEqual(self.case.knudsen_definition,
                         'lambda_over_source_distance')

    def test_metadata_records_kn_as_derived_and_disclaims_dsmc(self):
        provenance = self.results.provenance
        self.assertIn('derived metadata only',
                      provenance['knudsen']['role'])
        limitations = ' '.join(provenance['known_limitations'])
        self.assertIn('collisionless', limitations)
        self.assertIn('No DSMC data is read, written or compared',
                      limitations.replace('no DSMC', 'No DSMC'))

    # ---------------------------------------------------------- distributions
    def test_distribution_csv_covers_every_face_with_the_required_columns(self):
        path = self.case.surface_distribution_path
        self.assertTrue(path and os.path.isfile(path))
        with open(path, encoding='utf-8', newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), self.case.component_faces)
        for column in ('face_index', 'centroid_x', 'centroid_y', 'centroid_z',
                       'local_u', 'local_v', 'area', 'pressure',
                       'shear_stress', 'heat_flux', 'strike_count'):
            self.assertIn(column, rows[0])

    def test_distribution_coordinates_span_the_panel(self):
        with open(self.case.surface_distribution_path, encoding='utf-8',
                  newline='') as handle:
            rows = list(csv.DictReader(handle))
        local_u = np.array([float(row['local_u']) for row in rows])
        local_v = np.array([float(row['local_v']) for row in rows])
        # Every centroid lies inside the panel, and the centroids span it to
        # within one element (a triangle centroid sits inside its own cell,
        # so the extremes are inset by a fraction of the element size).
        self.assertLessEqual(float(np.max(np.abs(local_u))),
                             PANEL_LENGTH_U / 2)
        self.assertLessEqual(float(np.max(np.abs(local_v))),
                             PANEL_WIDTH_V / 2)
        element_u = PANEL_LENGTH_U / 44
        element_v = PANEL_WIDTH_V / 24
        self.assertGreater(float(np.ptp(local_u)),
                           PANEL_LENGTH_U - element_u)
        self.assertGreater(float(np.ptp(local_v)), PANEL_WIDTH_V - element_v)
        # Symmetric about the panel center.
        self.assertAlmostEqual(float(np.mean(local_u)), 0.0, places=6)
        self.assertAlmostEqual(float(np.mean(local_v)), 0.0, places=6)

    def test_integrated_pressure_force_matches_the_exported_faces(self):
        # Re-integrate the CSV independently: sum(-p_i * A_i * n_hat).
        with open(self.case.surface_distribution_path, encoding='utf-8',
                  newline='') as handle:
            rows = list(csv.DictReader(handle))
        pressure = np.array([float(row['pressure']) for row in rows])
        area = np.array([float(row['area']) for row in rows])
        recovered = -float(np.sum(pressure * area))       # along n
        np.testing.assert_allclose(self.case.pressure_force,
                                   [0.0, 0.0, recovered], rtol=1e-9,
                                   atol=1e-12)
        self.assertAlmostEqual(float(np.sum(area)), PANEL_AREA, places=6)

    def test_distribution_sidecar_records_units_and_provenance(self):
        sidecar = (os.path.splitext(self.case.surface_distribution_path)[0]
                   + '.meta.json')
        document = json.loads(Path(sidecar).read_text(encoding='utf-8'))
        self.assertEqual(document['units']['pressure'], 'Pa')
        self.assertEqual(document['units']['heat_flux'], 'W/m^2')
        self.assertEqual(document['plume_model'], 'SimplifiedGasKinetics')
        self.assertEqual(document['knudsen_number'], 0.25)
        self.assertEqual(document['source_axis_mode'], 'parallel_to_normal')
        self.assertIn('none', document['interpolation'])

    # --------------------------------------------------------------- artifacts
    def test_vtk_remains_the_primary_per_face_output(self):
        self.assertTrue(self.case.vtk_path
                        and os.path.isfile(self.case.vtk_path))
        self.assertTrue(self.case.vtk_path.endswith('.vtu'))

    def test_summary_csv_carries_the_new_columns(self):
        with open(self.results.summary_csv_path, encoding='utf-8',
                  newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        for column in ('source_offset_u', 'source_offset_v',
                       'source_axis_mode', 'model_variant', 'normal_force',
                       'local_moment_u', 'local_moment_v', 'local_moment_n',
                       'center_of_pressure_u', 'center_of_pressure_v',
                       'knudsen_number', 'mean_free_path',
                       'knudsen_reference_length', 'knudsen_definition',
                       'surface_distribution_path'):
            self.assertIn(column, row, f'missing CSV column {column}')
        self.assertEqual(row['model_variant'], 'Simplified')
        self.assertEqual(row['source_axis_mode'], 'parallel_to_normal')
        self.assertAlmostEqual(float(row['knudsen_number']), 0.25)

    def test_metadata_json_keeps_model_and_knudsen_provenance(self):
        document = json.loads(
            Path(self.results.metadata_path).read_text(encoding='utf-8'))
        provenance = document['provenance']
        self.assertEqual(provenance['plume_model'], 'SimplifiedGasKinetics')
        self.assertEqual(provenance['source_axis_mode'], 'parallel_to_normal')
        self.assertEqual(provenance['knudsen']['mean_free_path_m'], 1.0)
        basis = provenance['panel_basis']
        np.testing.assert_allclose(basis['u'], PANEL_U, atol=1e-12)
        np.testing.assert_allclose(basis['v'], PANEL_V, atol=1e-12)
        np.testing.assert_allclose(basis['n'], PANEL_N, atol=1e-12)
        self.assertEqual(document['cases'][0]['knudsen_number'], 0.25)

    def test_plot_generation_completes_headless(self):
        # setUpClass already ran plots via the configuration; re-running
        # returns the paths and must not raise without a display.
        written = self.study.plot()
        self.assertTrue(written)
        for path in written:
            self.assertTrue(os.path.isfile(path), path)
        names = {os.path.basename(path) for path in written}
        self.assertIn('normal_force_vs_offset_u.png', names)
        self.assertIn(f'panel_pressure_{self.case.case_id}.png', names)


class ModelSelectionChangesTheAnswer(unittest.TestCase):
    """The selected model drives the calculation, not just the metadata."""

    @classmethod
    def setUpClass(cls):
        cls.output_dir = tempfile.mkdtemp(prefix='pyrpod_iss_models_')
        cls.cases = {}
        for label, config_path in (('simplified', BASELINE_YAML),
                                   ('full_cai', FULL_CAI_YAML)):
            study = TradeStudy.from_config(
                config_path, output_dir=os.path.join(cls.output_dir, label))
            cls.cases[label] = study.run().cases[0]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.output_dir, ignore_errors=True)

    def test_both_models_run_the_same_geometry_and_pose(self):
        simplified, full = self.cases['simplified'], self.cases['full_cai']
        self.assertEqual(simplified.plume_model, 'SimplifiedGasKinetics')
        self.assertEqual(full.plume_model, 'CollisionlessGasKinetics')
        self.assertEqual(simplified.mesh_faces, full.mesh_faces)
        np.testing.assert_allclose(simplified.plume_source_position,
                                   full.plume_source_position, atol=1e-12)
        self.assertEqual(simplified.struck_faces, full.struck_faces)

    def test_the_loads_actually_differ(self):
        # If model selection were metadata only, these would be identical.
        simplified, full = self.cases['simplified'], self.cases['full_cai']
        self.assertNotAlmostEqual(simplified.max_pressure, full.max_pressure,
                                  places=6)
        self.assertNotAlmostEqual(simplified.normal_force, full.normal_force,
                                  places=6)

    def test_the_two_models_stay_physically_close(self):
        # Both are the same collisionless jet; the full model differs by the
        # near-field correction, not by an order of magnitude.
        simplified, full = self.cases['simplified'], self.cases['full_cai']
        for name in ('normal_force', 'max_pressure', 'max_heat_flux'):
            with self.subTest(quantity=name):
                a, b = getattr(simplified, name), getattr(full, name)
                self.assertLess(abs(a - b) / abs(a), 0.15)

    def test_both_report_the_same_derived_knudsen_number(self):
        # Kn depends on the configuration, never on the model.
        self.assertEqual(self.cases['simplified'].knudsen_number,
                         self.cases['full_cai'].knudsen_number)


class OffsetSweep(unittest.TestCase):
    """Offset-sweep enumeration and physics, on a small mesh."""

    DISTANCES = [3.0, 5.0]
    OFFSETS_U = [-6.0, 0.0, 6.0]

    @classmethod
    def setUpClass(cls):
        # A deliberately coarse panel: this test is about the sweep, not the
        # mesh. Same 22 x 12 m panel, 8 x 4 quads = 64 faces.
        import sys
        sys.path.insert(0, str(CASE_DIR / 'stl'))
        from generate_panel import build_panel_mesh

        cls.stl_path = CASE_DIR / 'stl' / 'iss_panel_coarse.stl'
        build_panel_mesh(n_u=8, n_v=4).save(str(cls.stl_path))

        cls.output_dir = tempfile.mkdtemp(prefix='pyrpod_iss_sweep_')
        mapping = coarse_panel_mapping(source_distances=cls.DISTANCES,
                                       source_offsets_u=cls.OFFSETS_U,
                                       source_offsets_v=[0.0])
        cls.study, cls.results = run_config(mapping, cls.output_dir, 'sweep')
        cls.cases = list(cls.results.cases)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.output_dir, ignore_errors=True)
        cls.stl_path.unlink(missing_ok=True)

    def _case(self, distance, offset_u):
        for case in self.cases:
            if (case.source_distance == distance
                    and case.source_offset_u == offset_u):
                return case
        raise AssertionError(f'no case at L={distance}, u={offset_u}')

    # ------------------------------------------------------------ enumeration
    def test_case_count_is_the_product_of_the_swept_axes(self):
        self.assertEqual(len(self.cases),
                         len(self.DISTANCES) * len(self.OFFSETS_U) * 1)

    def test_cases_are_ordered_distance_major_then_offset(self):
        self.assertEqual(
            [(case.source_distance, case.source_offset_u)
             for case in self.cases],
            [(distance, offset) for distance in self.DISTANCES
             for offset in self.OFFSETS_U])

    def test_case_ids_are_unique_and_name_their_parameters(self):
        ids = [case.case_id for case in self.cases]
        self.assertEqual(len(set(ids)), len(ids))
        self.assertEqual(self._case(3.0, -6.0).case_id,
                         'case000_modelSimplified_L3_um6_v0')
        self.assertEqual(self._case(5.0, 6.0).case_id,
                         'case005_modelSimplified_L5_u6_v0')

    def test_each_case_wrote_its_own_jfh_and_artifacts(self):
        for case in self.cases:
            with self.subTest(case=case.case_id):
                self.assertTrue(os.path.isfile(case.jfh_path))
                self.assertIn(case.case_id, case.jfh_path)
                self.assertTrue(os.path.isfile(case.vtk_path))
                self.assertTrue(
                    os.path.isfile(case.surface_distribution_path))

    # ---------------------------------------------------------------- physics
    def test_the_source_translates_without_tilting_the_plume_axis(self):
        for case in self.cases:
            with self.subTest(case=case.case_id):
                np.testing.assert_allclose(
                    case.plume_source_position,
                    [case.source_offset_u, case.source_offset_v,
                     case.source_distance], atol=1e-9)
                dcm = np.asarray(case.plume_source_orientation).reshape(3, 3)
                np.testing.assert_allclose(dcm[:, 0], -PANEL_N, atol=1e-9)

    def test_centered_cases_are_symmetric(self):
        for distance in self.DISTANCES:
            case = self._case(distance, 0.0)
            with self.subTest(L=distance):
                scale = case.normal_force * 11.0
                self.assertAlmostEqual(case.local_moment_v, 0.0,
                                       delta=1e-9 * scale)
                self.assertAlmostEqual(case.center_of_pressure_u, 0.0,
                                       places=6)

    def test_offset_source_moves_the_center_of_pressure_with_it(self):
        for distance in self.DISTANCES:
            for offset in (-6.0, 6.0):
                case = self._case(distance, offset)
                with self.subTest(L=distance, u=offset):
                    # Inside the panel and well clear of its edges, the
                    # footprint is fully captured, so the CoP sits at the
                    # centerline intersection.
                    self.assertAlmostEqual(case.center_of_pressure_u, offset,
                                           delta=0.15)
                    # Zero to within cross-product round-off on a 12 m span.
                    self.assertAlmostEqual(case.center_of_pressure_v, 0.0,
                                           delta=1e-3)

    def test_positive_u_offset_gives_a_positive_moment_about_v(self):
        # The documented sign convention (see pyrpod.mdao.surface_loads).
        for distance in self.DISTANCES:
            positive = self._case(distance, 6.0)
            negative = self._case(distance, -6.0)
            with self.subTest(L=distance):
                self.assertGreater(positive.local_moment_v, 0.0)
                self.assertLess(negative.local_moment_v, 0.0)
                # Mirror symmetry of the panel about u = 0.
                self.assertAlmostEqual(positive.local_moment_v,
                                       -negative.local_moment_v,
                                       delta=1e-6 * abs(positive.local_moment_v))
                self.assertAlmostEqual(positive.normal_force,
                                       negative.normal_force,
                                       delta=1e-6 * positive.normal_force)

    def test_peak_pressure_falls_with_stand_off(self):
        for offset in self.OFFSETS_U:
            with self.subTest(u=offset):
                self.assertGreater(self._case(3.0, offset).max_pressure,
                                   self._case(5.0, offset).max_pressure)

    def test_knudsen_tracks_the_swept_distance(self):
        for distance in self.DISTANCES:
            case = self._case(distance, 0.0)
            with self.subTest(L=distance):
                self.assertAlmostEqual(case.knudsen_number, 1.0 / distance,
                                       places=12)

    # ----------------------------------------------------------------- output
    def test_summary_csv_has_one_row_per_case_with_flat_columns(self):
        with open(self.results.summary_csv_path, encoding='utf-8',
                  newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), len(self.cases))
        offsets = sorted({float(row['source_offset_u']) for row in rows})
        self.assertEqual(offsets, sorted(self.OFFSETS_U))
        for row in rows:
            self.assertNotEqual(row['normal_force'], '')
            self.assertNotEqual(row['local_moment_v'], '')
            self.assertNotEqual(row['knudsen_number'], '')

    def test_plots_include_the_offset_trends(self):
        written = self.study.plot()
        names = {os.path.basename(path) for path in written}
        for expected in ('normal_force_vs_offset_u.png',
                         'moment_v_vs_offset_u.png',
                         'peak_pressure_vs_offset_u.png',
                         'cop_u_vs_offset_u.png',
                         'normal_force_vs_distance.png'):
            self.assertIn(expected, names)
        for path in written:
            self.assertTrue(os.path.isfile(path), path)


class TransverseOffsetSweep(unittest.TestCase):
    """A v offset loads the u moment, mirroring the u-offset behavior."""

    @classmethod
    def setUpClass(cls):
        import sys
        sys.path.insert(0, str(CASE_DIR / 'stl'))
        from generate_panel import build_panel_mesh

        cls.stl_path = CASE_DIR / 'stl' / 'iss_panel_coarse.stl'
        build_panel_mesh(n_u=8, n_v=4).save(str(cls.stl_path))

        cls.output_dir = tempfile.mkdtemp(prefix='pyrpod_iss_voffset_')
        mapping = coarse_panel_mapping(source_distances=[3.0],
                                       source_offsets_u=[0.0],
                                       source_offsets_v=[-3.0, 0.0, 3.0])
        cls.study, cls.results = run_config(mapping, cls.output_dir,
                                            'voffset')
        cls.cases = {case.source_offset_v: case for case in cls.results.cases}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.output_dir, ignore_errors=True)
        cls.stl_path.unlink(missing_ok=True)

    def test_three_cases_one_per_transverse_offset(self):
        self.assertEqual(sorted(self.cases), [-3.0, 0.0, 3.0])

    def test_positive_v_offset_gives_a_negative_moment_about_u(self):
        # M = (v_off * v_hat) x (-F_n * n_hat) = -v_off * F_n * u_hat.
        self.assertLess(self.cases[3.0].local_moment_u, 0.0)
        self.assertGreater(self.cases[-3.0].local_moment_u, 0.0)
        self.assertAlmostEqual(self.cases[0.0].local_moment_u, 0.0,
                               delta=1e-9 * self.cases[0.0].normal_force * 6.0)

    def test_center_of_pressure_follows_the_transverse_offset(self):
        for offset in (-3.0, 3.0):
            with self.subTest(v=offset):
                self.assertAlmostEqual(
                    self.cases[offset].center_of_pressure_v, offset,
                    delta=0.15)
                self.assertAlmostEqual(
                    self.cases[offset].center_of_pressure_u, 0.0, delta=1e-3)

    def test_transverse_plots_are_generated_when_v_is_swept(self):
        names = {os.path.basename(path) for path in self.study.plot()}
        for expected in ('normal_force_vs_offset_v.png',
                         'moment_u_vs_offset_v.png',
                         'cop_v_vs_offset_v.png'):
            self.assertIn(expected, names)


class BackwardCompatibility(unittest.TestCase):
    """Existing studies keep their behavior, paths and results."""

    def test_committed_flat_plate_configurations_still_parse(self):
        flat_plate = _TESTS_DIR.parent / 'case' / 'plume' / \
            'plume_flat_plate_sweep' / 'study'
        for name in ('flat_plate_baseline.yaml', 'flat_plate_sweep.yaml',
                     'flat_plate_sweep_single_jfh.yaml'):
            with self.subTest(config=name):
                config = StudyConfig.from_yaml(flat_plate / name)
                # Untouched defaults: the historical model, axis mode and
                # zero offsets, and no Knudsen metadata.
                self.assertEqual(config.plume_model, 'SimplifiedGasKinetics')
                self.assertEqual(config.sweep.source_axis_mode,
                                 'aim_at_reference')
                self.assertEqual(config.sweep.source_offsets_u, (0.0,))
                self.assertEqual(config.sweep.source_offsets_v, (0.0,))
                self.assertIsNone(config.knudsen)
                self.assertFalse(config.output.write_surface_distribution)

    def test_aim_at_reference_case_ids_are_unchanged(self):
        from pyrpod.mdao.plume_validation import case_id_for
        self.assertEqual(case_id_for(0, 0.0, 4.0), 'case000_alpha0p0_d4')
        self.assertEqual(case_id_for(12, -30.0, 2.0),
                         'case012_alpham30p0_d2')

    def test_baseline_flat_plate_results_are_unchanged(self):
        # The documented head-on baseline numbers of the existing case; the
        # new panel-local fields are additions, not a redefinition.
        flat_plate = (_TESTS_DIR.parent / 'case' / 'plume'
                      / 'plume_flat_plate_sweep' / 'study'
                      / 'flat_plate_baseline.yaml')
        output_dir = tempfile.mkdtemp(prefix='pyrpod_flat_plate_regression_')
        try:
            results = TradeStudy.from_config(
                flat_plate, output_dir=output_dir).run()
            case = results.cases[0]
            self.assertEqual(case.case_id, 'case000_alpha0p0_d4')
            np.testing.assert_allclose(case.force, [0.0, 0.0, -2.763],
                                       atol=5e-3)
            self.assertAlmostEqual(case.max_pressure, 0.340, places=3)
            self.assertAlmostEqual(case.max_heat_flux, 50.1, places=1)
            self.assertAlmostEqual(case.coefficients['CF'], 0.0391, places=4)
            # The additions are present and consistent with the old fields.
            self.assertAlmostEqual(case.normal_force, -float(case.force[2]),
                                   places=9)
            self.assertIsNone(case.knudsen_number)
            self.assertIsNone(case.surface_distribution_path)
        finally:
            shutil.rmtree(output_dir, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
