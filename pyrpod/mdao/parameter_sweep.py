"""
Single-history engine for prescribed plume/target parameter sweeps.

:class:`ParameterSweepStudy` runs a whole angle x distance sweep as ONE case
driven by ONE Jet Firing History: every pose contributes ``n_firings``
consecutive entries to a single history, the strike pipeline runs once over
all of them, and every strike file therefore belongs to the same run. It is
the ``sweep.mode: single_jfh`` engine behind the package-level API

>>> from pyrpod.mdao.TradeStudy import TradeStudy
>>> study = TradeStudy.from_config('sweep.yaml')     # mode: single_jfh
>>> results = study.run()

Why a second engine rather than a flag
--------------------------------------
The two decompositions answer different questions, and their outputs differ
in kind, not just in layout:

* :class:`pyrpod.mdao.plume_validation.PlumeValidationStudy` (per case) keeps
  poses independent -- the right thing when each pose is a separate
  experiment, as in a validation matrix compared pose-by-pose against
  reference data;
* this class treats the sweep as a single firing history -- the right thing
  when the poses form one sequence. The pipeline's cumulative fields then
  mean something: ``max_pressures`` / ``max_shears`` become the worst load
  each face saw ANYWHERE in the sweep, ``cum_strikes`` a coverage map, and
  ``cum_heat_flux_load`` the accumulated dose. Those envelopes are recorded
  in the study metadata and live in the VTK files, which land in one
  ``results/strikes/firing-<i>.vtu`` series that ParaView can scrub as a
  time sequence.

Everything else is shared with the per-case engine through
:mod:`pyrpod.mdao.study_runtime`: the same case objects, the same pose
generation, the same load integration, the same result schema, the same
reference comparison and plots. Per-firing numbers are identical between the
two engines -- only the decomposition into histories differs.

``n_firings`` semantics
-----------------------
``n_firings`` stays the number of entries EACH POSE contributes, so the same
configuration can be run either way. The single history is required to hold
exactly ``len(poses) * n_firings`` entries, checked after writing and again
after reading the file back.
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
from pyrpod.mdao.study_config import StudyConfig
from pyrpod.mdao.study_results import StudyResults
from pyrpod.mdao.study_runtime import TargetGeometry
from pyrpod.util.io.fs import ensure_dir

logger = logging.getLogger(__name__)

__all__ = ["ParameterSweepStudy"]


class ParameterSweepStudy:
    """Prescribed parameter sweep executed as one Jet Firing History."""

    def __init__(self, config: StudyConfig) -> None:
        self.config = config
        self.results: StudyResults | None = None
        self.comparison: ComparisonReport | None = None
        #: Worst-case per-face statistics over the whole sweep, per component
        #: (populated by run() when the full pipeline path is used).
        self.envelope: dict[str, dict[str, float]] = {}

    # ------------------------------------------------------------ lifecycle
    @classmethod
    def from_config(cls, path: str | os.PathLike[str]
                    ) -> "ParameterSweepStudy":
        """Build a sweep study from a YAML configuration file."""
        return cls(StudyConfig.from_yaml(path))

    @property
    def case_id(self) -> str:
        """The single case identifier every result row carries."""
        return self.config.study_name

    @property
    def jfh_name(self) -> str:
        """File name of the sweep's one Jet Firing History."""
        return f"{self.config.study_name}.A"

    # ------------------------------------------------------------------ run
    def run(self, write_outputs: bool = True) -> StudyResults:
        """Execute the whole sweep as one firing history.

        Returns
        -------
        StudyResults
            One :class:`~pyrpod.mdao.study_results.CaseResult` per firing and
            component, each keyed to the pose it realizes (``plate_angle_deg``
            and ``source_distance``) and sharing one ``case_id``.
        """
        config = self.config
        started = time.perf_counter()
        ensure_dir(config.output_dir)

        assets = study_runtime.load_case_assets(config)
        geometry = TargetGeometry.from_config(config, assets.target_vehicle.mesh)

        # 1. One history for the whole sweep, exactly poses x n_firings long.
        firings = firing_plan.build_sweep_firings(config.sweep, config.target)
        jfh_path = os.path.join(config.output_dir, "jfh", self.jfh_name)
        n_written = firing_plan.write_jfh_file(jfh_path, firings)
        expected = config.sweep.total_firings
        if n_written != expected:
            raise ValueError(
                f"sweep JFH {jfh_path!r} was written with {n_written} entries "
                f"but {len(config.sweep.poses)} poses x n_firings="
                f"{config.sweep.n_firings} requires exactly {expected}")

        jfh = study_runtime.read_generated_jfh(config, self.jfh_name)
        if len(jfh.JFH) != expected:
            raise ValueError(
                f"sweep JFH {jfh_path!r} reports {len(jfh.JFH)} firings, "
                f"expected {expected}")

        logger.info("Parameter sweep study started: study=%s mode=single_jfh "
                    "poses=%d firings_per_pose=%d jfh_entries=%d "
                    "components=%d target_faces=%d", config.study_name,
                    len(config.sweep.poses), config.sweep.n_firings, expected,
                    len(geometry.components), geometry.n_faces)

        # 2. One strike run over every firing in that history.
        firing_data, vtk_paths = study_runtime.compute_strikes(
            config, assets, geometry, jfh, output_root=config.output_dir)

        # 3. One result record per firing and component, keyed to its pose.
        results = StudyResults(
            study_name=config.study_name, output_dir=config.output_dir,
            provenance=study_runtime.study_provenance(
                config, geometry, jfh_path=os.path.abspath(jfh_path),
                case_id=self.case_id))
        timestamp = study_runtime.utc_timestamp()
        version = str(results.provenance.get("code_version", "unknown"))

        distribution_dir = os.path.join(
            config.output_dir, config.output.surface_distribution_subdir)

        for firing_index, firing in enumerate(firings):
            per_face = firing_data[str(firing_index + 1)]
            plate_angle_deg = _pose_value(firing.plate_angle_deg)
            source_distance = _pose_value(firing.source_distance)
            distributions = study_runtime.export_surface_distributions(
                config, geometry, firing, per_face, case_id=self.case_id,
                firing_id=firing_index + 1,
                plate_angle_deg=plate_angle_deg,
                source_distance=source_distance,
                output_dir=distribution_dir)
            results.cases.extend(study_runtime.build_case_results(
                config, geometry, firing, per_face,
                case_id=self.case_id, firing_id=firing_index + 1,
                plate_angle_deg=plate_angle_deg,
                source_distance=source_distance,
                jfh_path=jfh_path,
                vtk_path=vtk_paths[firing_index] if vtk_paths else None,
                code_version_id=version, timestamp=timestamp,
                surface_distribution_paths=distributions))

        # 4. The sweep envelope, which only a shared history can produce.
        last = firing_data[str(len(firings))]
        if all(name in last for name in study_runtime.CUMULATIVE_FIELDS):
            self.envelope = study_runtime.component_envelope(geometry, last)
            results.provenance["sweep_envelope"] = self.envelope
        else:
            self.envelope = {}
            logger.info("Sweep envelope unavailable: the cumulative fields "
                        "come from the full strike pipeline, which is "
                        "disabled when output.vtk.enabled is false.")

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

        logger.info("Parameter sweep study completed: study=%s firings=%d "
                    "records=%d wall_time_s=%.2f output=%s",
                    config.study_name, expected, len(results),
                    time.perf_counter() - started, config.output_dir)
        return results

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
        """Compare this sweep's results against a reference dataset.

        Records are matched on the swept parameters exactly as in the
        per-case engine, so the same reference data compares against either.
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
        """Generate the optional trend plots for this sweep."""
        results = results or self.results
        if results is None:
            raise RuntimeError("run() the study before plotting results")
        plot_dir = os.path.join(self.config.output_dir,
                                self.config.output.plots_subdir)
        return study_runtime.study_plots_for(self.config, results, plot_dir,
                                             comparison=self.comparison)


def _pose_value(value: float | None) -> float:
    """A firing's swept parameter; NaN when its pose was prescribed outright.

    Explicitly prescribed firings that carry no sweep parameterization are
    recorded as NaN rather than as a made-up angle or distance.
    """
    return float("nan") if value is None else float(value)
