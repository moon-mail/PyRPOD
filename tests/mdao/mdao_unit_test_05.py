# ========================
# PyRPOD: tests/mdao/mdao_unit_test_05.py
# ========================
# Unit tests for the YAML study configuration (pyrpod.mdao.study_config):
#
#   * the committed flat-plate example configurations parse, and every field
#     the study engine relies on comes back with the expected value;
#   * paths are resolved relative to the configuration file, and the case
#     must be a real PyRPOD case (config.ini present);
#   * validation rejects what must not be guessed: an unsupported plume
#     model, a bad firing count, a firing list that disagrees with
#     n_firings, non-orthogonal target axes, duplicate component names and
#     non-positive normalization values;
#   * incomplete normalization inputs disable the corresponding
#     coefficients instead of inventing defaults;
#   * backward compatibility: the study layer never touches the case's own
#     config.ini, which keeps parsing exactly as before.
#
# Run:  python -m pytest mdao/mdao_unit_test_05.py   (from tests/)

import configparser
import copy
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pytest
import yaml

from pyrpod.mdao.study_config import (
    Normalization,
    StudyConfig,
    StudyConfigError,
    SUPPORTED_PLUME_MODEL,
)

_TESTS_DIR = Path(__file__).resolve().parents[1]
CASE_DIR = _TESTS_DIR.parent / 'case' / 'plume' / 'plume_flat_plate_sweep'
BASELINE_YAML = CASE_DIR / 'study' / 'flat_plate_baseline.yaml'
SWEEP_YAML = CASE_DIR / 'study' / 'flat_plate_sweep.yaml'
SINGLE_JFH_YAML = CASE_DIR / 'study' / 'flat_plate_sweep_single_jfh.yaml'


def baseline_mapping():
    return yaml.safe_load(BASELINE_YAML.read_text(encoding='utf-8'))


def from_mapping(data):
    return StudyConfig.from_mapping(data, source_path=str(BASELINE_YAML))


class CommittedExampleConfigurations(unittest.TestCase):

    def test_baseline_configuration_parses(self):
        config = StudyConfig.from_yaml(BASELINE_YAML)

        self.assertEqual(config.study_name, 'cai2016_flat_plate_baseline')
        self.assertEqual(config.plume_model, SUPPORTED_PLUME_MODEL)
        self.assertTrue(os.path.isdir(config.case_dir))
        self.assertTrue(config.case_dir.endswith(os.sep))
        self.assertTrue(os.path.isfile(os.path.join(config.case_dir,
                                                    'config.ini')))
        self.assertEqual(config.source_path, str(BASELINE_YAML))

        self.assertEqual(config.sweep.plate_angles_deg, (0.0,))
        self.assertEqual(config.sweep.source_distances, (4.0,))
        self.assertEqual(config.sweep.n_firings, 1)
        self.assertEqual(config.sweep.thrusters, (1,))
        self.assertEqual(config.n_cases, 1)

        np.testing.assert_allclose(config.target.reference_point,
                                   [0.0, 0.0, 0.0])
        np.testing.assert_allclose(config.target.normal, [0.0, 0.0, 1.0])
        np.testing.assert_allclose(config.target.tangent, [1.0, 0.0, 0.0])
        self.assertEqual([c.name for c in config.target.components], ['plate'])

        np.testing.assert_allclose(config.loads.moment_reference_point,
                                   [0.0, 0.0, 0.0])
        self.assertTrue(config.loads.normalization.has_moment_inputs)
        self.assertTrue(config.output.write_vtk)
        self.assertFalse(config.output.write_plots)
        self.assertEqual(config.units['force'], 'N')

    def test_sweep_configuration_enumerates_every_combination(self):
        config = StudyConfig.from_yaml(SWEEP_YAML)
        self.assertEqual(len(config.sweep.plate_angles_deg), 19)
        self.assertEqual(len(config.sweep.source_distances), 5)
        self.assertEqual(config.n_cases, 95)
        self.assertTrue(config.output.write_plots)

    def test_provenance_records_the_configuration_source(self):
        config = StudyConfig.from_yaml(BASELINE_YAML)
        provenance = config.provenance()
        self.assertEqual(provenance['config_path'], str(BASELINE_YAML))
        self.assertEqual(provenance['plume_model'], SUPPORTED_PLUME_MODEL)
        self.assertIn('units', provenance)

    def test_output_directory_can_be_overridden(self):
        config = StudyConfig.from_yaml(BASELINE_YAML)
        with tempfile.TemporaryDirectory() as tmp:
            moved = config.with_output_dir(tmp)
            self.assertEqual(moved.output_dir, os.path.abspath(tmp))
            # The original is untouched (the dataclass is frozen).
            self.assertNotEqual(config.output_dir, moved.output_dir)


class PathResolution(unittest.TestCase):

    def test_relative_paths_resolve_against_the_configuration_file(self):
        config = StudyConfig.from_yaml(BASELINE_YAML)
        self.assertEqual(Path(config.case_dir).resolve(), CASE_DIR.resolve())

    def test_missing_case_directory_is_rejected(self):
        data = baseline_mapping()
        data['study']['case_dir'] = 'no/such/case'
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('case_dir', str(excinfo.value))

    def test_directory_without_config_ini_is_rejected(self):
        data = baseline_mapping()
        with tempfile.TemporaryDirectory() as tmp:
            data['study']['case_dir'] = tmp
            with pytest.raises(StudyConfigError) as excinfo:
                from_mapping(data)
            self.assertIn('config.ini', str(excinfo.value))

    def test_missing_configuration_file_is_reported(self):
        with pytest.raises(StudyConfigError):
            StudyConfig.from_yaml(CASE_DIR / 'study' / 'does_not_exist.yaml')


class Validation(unittest.TestCase):

    def test_study_name_is_required(self):
        data = baseline_mapping()
        del data['study']['name']
        with pytest.raises(StudyConfigError):
            from_mapping(data)

    def test_every_supported_plume_model_is_accepted(self):
        # Both collisionless Cai variants may be selected by name; the model
        # named here is the one that computes the plume field.
        for name in ('SimplifiedGasKinetics', 'CollisionlessGasKinetics'):
            with self.subTest(model=name):
                data = baseline_mapping()
                data['plume_model']['name'] = name
                self.assertEqual(from_mapping(data).plume_model, name)

    def test_unknown_plume_model_is_rejected(self):
        data = baseline_mapping()
        data['plume_model']['name'] = 'DSMCGasKinetics'
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('DSMCGasKinetics', str(excinfo.value))
        self.assertIn(SUPPORTED_PLUME_MODEL, str(excinfo.value))

    def test_omitted_plume_model_keeps_the_historical_default(self):
        data = baseline_mapping()
        data.pop('plume_model')
        self.assertEqual(from_mapping(data).plume_model, SUPPORTED_PLUME_MODEL)

    def test_invalid_firing_counts_are_rejected(self):
        for value in (0, -2, 1.5, 'many'):
            with self.subTest(value=value):
                data = baseline_mapping()
                data['sweep']['n_firings'] = value
                with pytest.raises(StudyConfigError):
                    from_mapping(data)

    def test_explicit_firings_must_match_n_firings(self):
        data = baseline_mapping()
        data['sweep']['n_firings'] = 3
        data['sweep']['firings'] = [
            {'position': [0.0, 0.0, 4.0], 'dcm': np.eye(3).tolist()}]
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('n_firings', str(excinfo.value))

    def test_explicit_firings_are_parsed_when_the_count_agrees(self):
        data = baseline_mapping()
        data['sweep']['n_firings'] = 2
        data['sweep']['firings'] = [
            {'position': [0.0, 0.0, 4.0], 'dcm': np.eye(3).tolist(),
             'thrusters': [1], 'duration_s': 0.25},
            {'position': [0.5, 0.0, 4.0], 'dcm': np.eye(3).tolist()}]
        config = from_mapping(data)
        self.assertEqual(len(config.sweep.firings), 2)
        self.assertEqual(config.sweep.firings[0].duration_s, 0.25)
        # The second firing inherits the sweep defaults.
        self.assertEqual(config.sweep.firings[1].duration_s, 1.0)
        self.assertEqual(config.sweep.firings[1].thrusters, (1,))

    def test_malformed_firing_pose_is_rejected(self):
        data = baseline_mapping()
        data['sweep']['n_firings'] = 1
        data['sweep']['firings'] = [
            {'position': [0.0, 0.0, 4.0], 'dcm': [[1.0, 0.0], [0.0, 1.0]]}]
        with pytest.raises(StudyConfigError):
            from_mapping(data)

    def test_target_axes_must_be_orthogonal(self):
        data = baseline_mapping()
        data['target']['tangent'] = [0.0, 0.0, 1.0]
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('orthogonal', str(excinfo.value))

    def test_empty_sweep_axes_are_rejected(self):
        for key in ('plate_angles_deg', 'source_distances'):
            with self.subTest(key=key):
                data = baseline_mapping()
                data['sweep'][key] = []
                with pytest.raises(StudyConfigError):
                    from_mapping(data)

    def test_non_positive_distances_are_rejected(self):
        data = baseline_mapping()
        data['sweep']['source_distances'] = [4.0, 0.0]
        with pytest.raises(StudyConfigError):
            from_mapping(data)

    def test_duplicate_component_names_are_rejected(self):
        data = baseline_mapping()
        data['target']['components'] = [{'name': 'plate'}, {'name': 'plate'}]
        with pytest.raises(StudyConfigError):
            from_mapping(data)

    def test_component_bounds_are_validated(self):
        data = baseline_mapping()
        data['target']['components'] = [
            {'name': 'plate',
             'bounds': {'min': [1.0, 1.0, 1.0], 'max': [0.0, 0.0, 0.0]}}]
        with pytest.raises(StudyConfigError):
            from_mapping(data)


class SweepMode(unittest.TestCase):
    """per_case (default) vs single_jfh, and what n_firings means in each."""

    def test_default_mode_is_per_case(self):
        self.assertEqual(StudyConfig.from_yaml(BASELINE_YAML).sweep.mode,
                         'per_case')

    def test_single_jfh_configuration_parses(self):
        config = StudyConfig.from_yaml(SINGLE_JFH_YAML)
        self.assertEqual(config.sweep.mode, 'single_jfh')
        self.assertEqual(config.n_cases, 95)
        self.assertEqual(config.sweep.n_firings, 1)
        # n_firings is per pose; the one history holds poses x n_firings.
        self.assertEqual(config.sweep.total_firings, 95)

    def test_unknown_mode_is_rejected(self):
        data = baseline_mapping()
        data['sweep']['mode'] = 'one_big_jfh'
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('sweep.mode', str(excinfo.value))

    def test_total_firings_scales_with_firings_per_pose(self):
        data = baseline_mapping()
        data['sweep']['mode'] = 'single_jfh'
        data['sweep']['plate_angles_deg'] = [-30.0, 0.0, 30.0]
        data['sweep']['source_distances'] = [2.0, 4.0]
        data['sweep']['n_firings'] = 4
        sweep = from_mapping(data).sweep
        self.assertEqual(len(sweep.poses), 6)
        self.assertEqual(sweep.total_firings, 24)

    def test_poses_are_enumerated_distance_major(self):
        data = baseline_mapping()
        data['sweep']['plate_angles_deg'] = [-10.0, 10.0]
        data['sweep']['source_distances'] = [2.0, 4.0]
        self.assertEqual(from_mapping(data).sweep.poses,
                         ((-10.0, 2.0), (10.0, 2.0),
                          (-10.0, 4.0), (10.0, 4.0)))

    def test_explicit_firing_count_is_checked_per_mode(self):
        one_firing = [{'position': [0.0, 0.0, 4.0],
                       'dcm': np.eye(3).tolist()}]

        # per_case: one firing per case is enough for any number of poses.
        data = baseline_mapping()
        data['sweep']['plate_angles_deg'] = [-10.0, 10.0]
        data['sweep']['n_firings'] = 1
        data['sweep']['firings'] = one_firing
        self.assertEqual(len(from_mapping(data).sweep.firings), 1)

        # single_jfh: the sequence must cover every pose.
        data['sweep']['mode'] = 'single_jfh'
        with pytest.raises(StudyConfigError) as excinfo:
            from_mapping(data)
        self.assertIn('single_jfh', str(excinfo.value))

        data['sweep']['firings'] = one_firing * 2
        self.assertEqual(len(from_mapping(data).sweep.firings), 2)


class CoefficientNormalization(unittest.TestCase):

    def test_complete_inputs_enable_force_and_moment_coefficients(self):
        normalization = StudyConfig.from_yaml(
            BASELINE_YAML).loads.normalization
        self.assertTrue(normalization.has_force_inputs)
        self.assertTrue(normalization.has_moment_inputs)
        self.assertEqual(normalization.reference_area, 64.0)
        self.assertEqual(normalization.reference_length, 4.0)

    def test_absent_normalization_leaves_every_input_unset(self):
        data = baseline_mapping()
        del data['loads']['normalization']
        normalization = from_mapping(data).loads.normalization
        self.assertFalse(normalization.has_force_inputs)
        self.assertFalse(normalization.has_moment_inputs)
        self.assertEqual(normalization.to_dict(),
                         {'reference_area': None, 'reference_length': None,
                          'dynamic_pressure': None,
                          'reference_heat_flux': None})

    def test_partial_normalization_disables_only_what_it_must(self):
        data = baseline_mapping()
        del data['loads']['normalization']['reference_length']
        normalization = from_mapping(data).loads.normalization
        self.assertTrue(normalization.has_force_inputs)
        self.assertFalse(normalization.has_moment_inputs)

    def test_non_positive_values_are_rejected(self):
        with pytest.raises(StudyConfigError):
            Normalization.from_mapping({'dynamic_pressure': 0.0})


class BackwardCompatibility(unittest.TestCase):
    """The study layer sits on top of the case; it changes nothing below."""

    def test_case_config_ini_still_parses_unchanged(self):
        config = configparser.ConfigParser()
        config.read(str(CASE_DIR / 'config.ini'))
        self.assertEqual(config['pm']['kinetics'], 'Simplified')
        self.assertEqual(config['tv']['stl'], 'flat_plate_transformed.stl')
        self.assertEqual(config['jfh']['jfh'], 'jfh_flat_plate_sweep.A')

    def test_geometry_id_defaults_to_the_case_target_stl(self):
        data = baseline_mapping()
        del data['target']['geometry_id']
        self.assertEqual(from_mapping(data).target.geometry_id,
                         'flat_plate_transformed.stl')

    def test_defaults_apply_to_a_minimal_configuration(self):
        minimal = {'study': {'name': 'minimal',
                             'case_dir': str(CASE_DIR)},
                   'sweep': {'plate_angles_deg': [0.0],
                             'source_distances': [4.0]}}
        config = StudyConfig.from_mapping(copy.deepcopy(minimal),
                                          source_path=str(BASELINE_YAML))
        self.assertEqual(config.sweep.n_firings, 1)
        self.assertEqual(config.plume_model, SUPPORTED_PLUME_MODEL)
        self.assertEqual([c.name for c in config.target.components],
                         ['target'])
        self.assertFalse(config.loads.normalization.has_force_inputs)
        self.assertTrue(config.output.write_vtk)


if __name__ == '__main__':
    unittest.main()
