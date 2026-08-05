"""
Per-case engine for prescribed plume/target validation studies.

:class:`PlumeValidationStudy` runs a configured sweep with ONE Jet Firing
History per angle-distance combination: every pose is an independent case
with its own history, its own strike run and its own artifacts, so no pose's
results can be polluted by another's. It is the ``sweep.mode: per_case``
(default) engine behind the package-level API

>>> from pyrpod.mdao.TradeStudy import TradeStudy
>>> study = TradeStudy.from_config('case.yaml')
>>> results = study.run()

The sibling :class:`pyrpod.mdao.parameter_sweep.ParameterSweepStudy` runs the
same configuration as a SINGLE history spanning every pose
(``sweep.mode: single_jfh``); both share the plumbing in
:mod:`pyrpod.mdao.study_runtime`, so their per-firing numbers are identical
and only the decomposition into histories differs.

Per case this engine does exactly seven things: build the prescribed firings,
write the case's Jet Firing History, run the plume-strike calculation,
integrate per-component loads, export the per-face VTK fields, record a
structured result, and (optionally) compare against external reference data.

Existing PyRPOD objects do all the domain work -- ``JetFiringHistory``,
``TargetVehicle``, ``VisitingVehicle``, ``MissionEnvironment`` and
``PlumeStrikeEstimationStudy`` -- so a study inherits the pipeline's
validation, logging, plume physics and VTK conventions rather than
reimplementing them.

Known limitation
----------------
Face selection is the pipeline's existing one: a face is struck when it lies
inside the plume wedge/radius and its normal faces the source. Plume
shadowing, occlusion by intervening geometry and self-shadowing of concave
targets are NOT modeled, so a face hidden behind other geometry still
receives load. This branch does not change that behavior; see
``docs/plume_validation_study.md``.
"""

from __future__ import annotations

import logging
import os
import time

from pyrpod.mdao import firing_plan, study_runtime
from pyrpod.mdao.reference_data import (
    ComparisonReport,
    ReferenceDataset,
    compare_results,
    load_reference_dataset,
)
from pyrpod.mdao.study_config import StudyConfig, SweepPose
from pyrpod.mdao.study_results import StudyResults
from pyrpod.mdao.study_runtime import (
    CaseAssets,
    TargetGeometry,
    study_plots_for,
)
from pyrpod.util.io.fs import ensure_dir

logger = logging.getLogger(__name__)

__all__ = ["PlumeValidationStudy", "case_id_for"]


def _token(value: float, spec: str = ".4g") -> str:
    """Filesystem-safe number token: '-' -> 'm', '.' -> 'p'."""
    return format(value, spec).replace("-", "m").replace(".", "p")


def case_id_for(index: int, plate_angle_deg: float, source_distance: float,
                pose: SweepPose | None = None,
                model_variant: str | None = None) -> str:
    """Stable, filesystem-safe case identifier.

    The identifier names the parameters that actually vary in the study's
    axis mode, so a case directory is self-describing:

    * ``aim_at_reference`` (default) keeps the historical
      ``case000_alpha0p0_d4`` form, byte-for-byte, so existing studies write
      to exactly the paths they always did;
    * ``parallel_to_normal`` names the model, stand-off and panel-local
      offsets instead, e.g. ``case000_modelCollisionless_L4_u0_v0`` -- the
      approach angle is fixed in that mode and would carry no information.
    """
    if pose is not None and pose.axis_mode == "parallel_to_normal":
        model = f"_model{model_variant}" if model_variant else ""
        return (f"case{index:03d}{model}"
                f"_L{_token(pose.source_distance)}"
                f"_u{_token(pose.source_offset_u)}"
                f"_v{_token(pose.source_offset_v)}")
    return (f"case{index:03d}_alpha{_token(plate_angle_deg, '.1f')}"
            f"_d{_token(source_distance)}")


class PlumeValidationStudy:
    """Prescribed validation sweep with one Jet Firing History per case."""

    def __init__(self, config: StudyConfig) -> None:
        self.config = config
        self.results: StudyResults | None = None
        self.comparison: ComparisonReport | None = None

    # ------------------------------------------------------------ lifecycle
    @classmethod
    def from_config(cls, path: str | os.PathLike[str]
                    ) -> "PlumeValidationStudy":
        """Build a study from a YAML configuration file."""
        return cls(StudyConfig.from_yaml(path))

    # ------------------------------------------------------------------ run
    def run(self, write_outputs: bool = True) -> StudyResults:
        """Execute every angle x distance case and return structured results.

        Parameters
        ----------
        write_outputs : bool, optional
            When True (default) the summary CSV, the JSON metadata document
            and the optional plots are written to the study output directory.
            The per-case VTK export is controlled by the configuration's own
            ``output.vtk.enabled`` flag.

        Returns
        -------
        StudyResults
            One :class:`~pyrpod.mdao.study_results.CaseResult` per case,
            component and firing.
        """
        config = self.config
        started = time.perf_counter()
        ensure_dir(config.output_dir)

        assets = study_runtime.load_case_assets(config)
        geometry = TargetGeometry.from_config(config, assets.target_vehicle.mesh)

        results = StudyResults(
            study_name=config.study_name, output_dir=config.output_dir,
            provenance=study_runtime.study_provenance(config, geometry))

        logger.info("Plume validation study started: study=%s mode=per_case "
                    "cases=%d firings_per_case=%d components=%d "
                    "target_faces=%d", config.study_name, config.n_cases,
                    config.sweep.n_firings, len(geometry.components),
                    geometry.n_faces)

        model_variant = config.plume_model.replace("GasKinetics", "")
        for index, pose in enumerate(config.sweep.sweep_poses):
            self._run_case(
                case_id=case_id_for(index, pose.plate_angle_deg,
                                    pose.source_distance, pose=pose,
                                    model_variant=model_variant),
                pose_index=index, pose=pose, results=results,
                assets=assets, geometry=geometry)

        if write_outputs:
            self.write_outputs(results)

        self.results = results
        if config.reference.path:
            self.comparison = self.compare(
                load_reference_dataset(config.reference.path,
                                       label=config.reference.label))
            if write_outputs and len(self.comparison):
                self.comparison.write_csv(os.path.join(
                    config.output_dir, "reference_comparison.csv"))

        logger.info("Plume validation study completed: study=%s cases=%d "
                    "records=%d wall_time_s=%.2f output=%s",
                    config.study_name, config.n_cases, len(results),
                    time.perf_counter() - started, config.output_dir)
        return results

    # ------------------------------------------------------------ one case
    def _run_case(self, *, case_id: str, pose_index: int, pose: SweepPose,
                  results: StudyResults, assets: CaseAssets,
                  geometry: TargetGeometry) -> None:
        config = self.config
        plate_angle_deg = pose.plate_angle_deg
        source_distance = pose.source_distance
        firings = firing_plan.build_case_firings(
            config.sweep, config.target, plate_angle_deg, source_distance,
            pose_index=pose_index, pose=pose)

        jfh_path = os.path.join(config.output_dir, "jfh", f"{case_id}.A")
        n_written = firing_plan.write_jfh_file(jfh_path, firings)
        if n_written != config.sweep.n_firings:
            raise ValueError(
                f"case {case_id}: wrote {n_written} JFH entries for a "
                f"requested n_firings of {config.sweep.n_firings}")

        jfh = study_runtime.read_generated_jfh(config, f"{case_id}.A")
        if len(jfh.JFH) != config.sweep.n_firings:
            raise ValueError(
                f"case {case_id}: JFH file {jfh_path!r} reports "
                f"{len(jfh.JFH)} firings, expected {config.sweep.n_firings}")

        firing_data, vtk_paths = study_runtime.compute_strikes(
            config, assets, geometry, jfh,
            output_root=os.path.join(config.output_dir, "cases", case_id))

        timestamp = study_runtime.utc_timestamp()
        version = str(results.provenance.get("code_version", "unknown"))

        distribution_dir = os.path.join(
            config.output_dir, "cases", case_id,
            config.output.surface_distribution_subdir)

        for firing_index, firing in enumerate(firings):
            per_face = firing_data[str(firing_index + 1)]
            distributions = study_runtime.export_surface_distributions(
                config, geometry, firing, per_face, case_id=case_id,
                firing_id=firing_index + 1, plate_angle_deg=plate_angle_deg,
                source_distance=source_distance,
                output_dir=distribution_dir)
            results.cases.extend(study_runtime.build_case_results(
                config, geometry, firing, per_face,
                case_id=case_id, firing_id=firing_index + 1,
                plate_angle_deg=plate_angle_deg,
                source_distance=source_distance, jfh_path=jfh_path,
                vtk_path=vtk_paths[firing_index] if vtk_paths else None,
                code_version_id=version, timestamp=timestamp,
                surface_distribution_paths=distributions))

    # -------------------------------------------------------------- outputs
    def write_outputs(self, results: StudyResults) -> None:
        """Write the summary CSV, the JSON metadata and any configured plots."""
        config = self.config
        results.write_csv(os.path.join(config.output_dir,
                                       config.output.summary_csv))
        results.write_metadata(os.path.join(config.output_dir,
                                            config.output.summary_metadata))
        if config.output.write_plots:
            results.plot_paths = self.plot(results)

    # ----------------------------------------------------------- comparison
    def compare(self, dataset: ReferenceDataset) -> ComparisonReport:
        """Compare this study's results against a reference dataset.

        The dataset is opaque: DSMC, analytical, experimental or another
        code, all handled identically (see :mod:`pyrpod.mdao.reference_data`).
        """
        if self.results is None:
            raise RuntimeError("run() the study before comparing results")
        report = compare_results(
            self.results.cases, dataset,
            reference_length=self.config.loads.normalization.reference_length)
        self.comparison = report
        return report

    # ---------------------------------------------------------------- plots
    def plot(self, results: StudyResults | None = None) -> list[str]:
        """Generate the optional trend plots for this study.

        Plot generation is entirely optional -- automated tests never require
        it -- and the plotting modules are imported lazily so a headless run
        pays no matplotlib cost unless plots were asked for.
        """
        results = results or self.results
        if results is None:
            raise RuntimeError("run() the study before plotting results")
        plot_dir = os.path.join(self.config.output_dir,
                                self.config.output.plots_subdir)
        return study_plots_for(self.config, results, plot_dir,
                               comparison=self.comparison)
