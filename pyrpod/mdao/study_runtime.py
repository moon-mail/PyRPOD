"""
Shared runtime for prescribed plume/target validation studies.

Both study styles -- :class:`pyrpod.mdao.plume_validation.PlumeValidationStudy`
(one Jet Firing History per angle-distance case) and
:class:`pyrpod.mdao.parameter_sweep.ParameterSweepStudy` (one history spanning
the whole sweep) -- need exactly the same plumbing: build the case objects,
precompute the target geometry, read back a generated JFH, run the strike
calculation, and turn one firing's per-face arrays into result records. That
plumbing lives here so the two classes stay siblings with no duplicated logic
and no inheritance between them.

Nothing in this module decides HOW a sweep is decomposed into histories; that
is the studies' own business.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import numpy as np
from numpy.typing import NDArray

from pyrpod.mdao.firing_plan import Firing
from pyrpod.mdao.study_config import StudyConfig
from pyrpod.mdao.study_results import CaseResult, code_version
from pyrpod.mdao.surface_distribution import (
    distribution_rows,
    write_surface_distribution,
)
from pyrpod.mdao.surface_loads import (
    face_areas,
    flow_directions,
    integrate_component_loads,
    project_to_panel_frame,
    select_component_faces,
)
from pyrpod.mission import MissionEnvironment
from pyrpod.plume.PlumeStrikeCalculator import (
    compute_face_centroids,
    compute_plume_strikes,
)
from pyrpod.plume.gas_kinetics_models import KINETICS_DISABLED, kinetics_key_for
from pyrpod.rpod import JetFiringHistory, PlumeStrikeEstimationStudy
from pyrpod.util.io.fs import ensure_dir
from pyrpod.vehicle import TargetVehicle, VisitingVehicle

logger = logging.getLogger(__name__)

__all__ = [
    "CaseAssets",
    "TargetGeometry",
    "build_case_results",
    "component_envelope",
    "compute_strikes",
    "export_surface_distributions",
    "knudsen_metadata",
    "load_case_assets",
    "read_generated_jfh",
    "study_plots_for",
    "study_provenance",
    "utc_timestamp",
]

#: Cumulative per-face fields the strike pipeline accumulates across a single
#: Jet Firing History. They are a sweep ENVELOPE only when one history spans
#: the sweep; with one history per case they simply restate that case.
CUMULATIVE_FIELDS = ("cum_strikes", "max_pressures", "max_shears",
                     "cum_heat_flux_load")


@dataclass
class CaseAssets:
    """The existing PyRPOD objects a study runs on."""

    target_vehicle: Any
    visiting_vehicle: Any
    environment: Any


@dataclass
class TargetGeometry:
    """Per-face target geometry, computed once and reused by every firing."""

    mesh: Any
    centroids: NDArray[np.float64]
    normals: NDArray[np.float64]
    areas: NDArray[np.float64]
    n_faces: int
    components: list[tuple[str, NDArray[np.int64]]]

    @classmethod
    def from_config(cls, config: StudyConfig, mesh: Any) -> "TargetGeometry":
        centroids = compute_face_centroids(mesh.vectors)
        n_faces = int(len(mesh.vectors))
        components = [
            (component.name,
             select_component_faces(component, centroids, n_faces))
            for component in config.target.components
        ]
        return cls(mesh=mesh, centroids=centroids,
                   normals=mesh.get_unit_normals(),
                   areas=face_areas(mesh.vectors), n_faces=n_faces,
                   components=components)


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_case_assets(config: StudyConfig) -> CaseAssets:
    """Construct the existing PyRPOD case objects for a study.

    Fails early and specifically when the case does not support the study
    workflow: gas kinetics disabled, or a configured thruster id that the
    case's thruster configuration file does not define.

    The study's ``plume_model`` selects which collisionless model computes
    the plume field. It is applied by setting the environment's IN-MEMORY
    ``[pm] kinetics`` key, which is the one input both strike paths read, so
    the selection reaches the calculation itself rather than only the
    metadata. The case's ``config.ini`` on disk is never modified, and a
    study naming the model its case already configures changes nothing.
    """
    case_dir = config.case_dir

    target_vehicle = TargetVehicle.TargetVehicle(case_dir)
    target_vehicle.set_stl()

    visiting_vehicle = VisitingVehicle.VisitingVehicle(case_dir)
    visiting_vehicle.set_thruster_config()
    visiting_vehicle.set_thruster_metrics()

    thruster_id = config.thruster_id
    if thruster_id is not None:
        available = list(getattr(visiting_vehicle, "thruster_data", {}) or {})
        if thruster_id not in available:
            raise ValueError(
                f"thruster.id {thruster_id!r} is not configured in the case's "
                f"thruster configuration file; available: {available}")

    environment = MissionEnvironment.MissionEnvironment(case_dir)
    configured = environment.config["pm"]["kinetics"]
    if configured == KINETICS_DISABLED:
        raise ValueError(
            f"case {case_dir!r} sets [pm] kinetics = {KINETICS_DISABLED}, "
            "which disables the gas-kinetics surface loads; a validation "
            "study needs pressure, shear and heat flux, so configure "
            "'Simplified' or 'Collisionless'")

    selected = kinetics_key_for(config.plume_model)
    if selected != configured:
        logger.info("Study selects plume model %s: [pm] kinetics %r -> %r "
                    "for this run (case config.ini is not modified)",
                    config.plume_model, configured, selected)
        environment.config["pm"]["kinetics"] = selected

    return CaseAssets(target_vehicle=target_vehicle,
                      visiting_vehicle=visiting_vehicle,
                      environment=environment)


def study_provenance(config: StudyConfig, geometry: TargetGeometry,
                     **extra: Any) -> dict[str, Any]:
    """Study-level provenance recorded in the metadata document."""
    u_hat, v_hat, n_hat = config.target.local_basis()
    provenance = config.provenance()
    provenance.update({
        "code_version": code_version(
            os.path.dirname(os.path.abspath(__file__))),
        "mesh_faces": geometry.n_faces,
        "geometry_id": config.target.geometry_id,
        "sweep_mode": config.sweep.mode,
        "n_cases": config.n_cases,
        "n_firings_per_pose": config.sweep.n_firings,
        "total_firings": config.sweep.total_firings,
        "components": [name for name, _ in geometry.components],
        "panel_basis": {
            "u": [float(value) for value in u_hat],
            "v": [float(value) for value in v_hat],
            "n": [float(value) for value in n_hat],
            "origin": [float(value)
                       for value in config.target.reference_point],
            "convention": "u = target tangent (longitudinal), "
                          "v = n x u (transverse), n = target normal "
                          "(toward the plume source); u x v = n. "
                          "normal_force = -F.n is positive INTO the panel.",
        },
        "known_limitations": [
            "plume shadowing, occlusion and back-facing geometry are not "
            "modeled (existing pipeline face-selection behavior)",
            "the plume models are collisionless: any configured Knudsen "
            "number is derived metadata and never enters the solution",
            "no DSMC data is read, written or compared here",
        ],
    })
    provenance.update(extra)
    return provenance


def read_generated_jfh(config: StudyConfig, jfh_name: str) -> Any:
    """Load a generated JFH through the normal JetFiringHistory reader.

    The JFH object keeps the case's ``config.ini`` (so every other setting is
    the case's own) but resolves its JFH asset from the study's own output
    tree, which is where generated firing histories live. ``case_dir`` is the
    documented asset-resolution root of ``JetFiringHistory``, so pointing it
    at the study output directory is exactly its intended use.
    """
    jfh = JetFiringHistory.JetFiringHistory(config.case_dir)
    jfh.case_dir = config.output_dir.rstrip("\\/") + os.sep
    if not jfh.config.has_section("jfh"):
        jfh.config.add_section("jfh")
    jfh.config.set("jfh", "jfh", jfh_name)
    jfh.read_jfh()
    if getattr(jfh, "JFH", None) is None:
        raise ValueError(
            f"failed to read generated JFH {jfh_name!r} from "
            f"{os.path.join(config.output_dir, 'jfh')}")
    return jfh


def compute_strikes(config: StudyConfig, assets: CaseAssets,
                    geometry: TargetGeometry, jfh: Any,
                    output_root: str) -> tuple[dict[str, dict[str, Any]],
                                               list[str]]:
    """Run one Jet Firing History through the plume-strike calculation.

    With VTK output enabled the full
    ``PlumeStrikeEstimationStudy.jfh_plume_strikes()`` pipeline runs, writing
    the standard per-face strike files to
    ``<output_root>/results/strikes/firing-<i>.vtu``; the target vehicle's
    output root is redirected there for the duration of the run so studies
    that execute several histories cannot overwrite one another's artifacts.
    With VTK disabled the pipeline's own per-firing core runs instead --
    identical numbers, no files.

    Returns
    -------
    (dict, list of str)
        The per-firing field dictionary keyed ``'1'``..``'N'``, and the
        per-firing VTK paths (empty when VTK output is disabled).
    """
    if not config.output.write_vtk:
        firing_data: dict[str, dict[str, Any]] = {}
        for index in range(len(jfh.JFH)):
            step = {"thrusters": jfh.JFH[index]["thrusters"],
                    "xyz": np.array(jfh.JFH[index]["xyz"]),
                    "dcm": np.array(jfh.JFH[index]["dcm"]),
                    "t": float(jfh.JFH[index]["t"])}
            result = compute_plume_strikes(
                geometry.mesh, geometry.normals, assets.visiting_vehicle,
                step, assets.environment, face_centroids=geometry.centroids)
            firing_data[str(index + 1)] = dict(result)
        return firing_data, []

    study = PlumeStrikeEstimationStudy.PlumeStrikeEstimationStudy(
        assets.environment)
    study.study_init(jfh, assets.target_vehicle, assets.visiting_vehicle)

    target_vehicle = assets.target_vehicle
    original_case_dir = target_vehicle.case_dir
    redirected = output_root.rstrip("\\/") + os.sep
    ensure_dir(os.path.join(redirected, "results", "strikes"))
    try:
        target_vehicle.case_dir = redirected
        firing_data = study.jfh_plume_strikes()
    finally:
        target_vehicle.case_dir = original_case_dir

    vtk_paths = [os.path.join(redirected, "results", "strikes",
                              f"firing-{index}.vtu")
                 for index in range(len(jfh.JFH))]
    return firing_data, vtk_paths


def knudsen_metadata(config: StudyConfig,
                     source_distance: float) -> dict[str, Any]:
    """Derived Knudsen fields for one case, or empty when unconfigured.

    Kn is METADATA: it is computed from the mean free path the configuration
    supplied and the reference length it selected, and no part of the
    analytical plume solution reads it back. When the study has no
    ``knudsen`` block every Knudsen field is simply left unset.
    """
    spec = config.knudsen
    if spec is None:
        return {}
    metadata: dict[str, Any] = {
        "mean_free_path": spec.mean_free_path_m,
        "knudsen_definition": spec.definition,
    }
    if (spec.reference_mode == "source_distance"
            and not np.isfinite(source_distance)):
        # An explicitly prescribed firing carries no swept distance (recorded
        # as NaN), so a distance-referenced Kn does not exist for it. The
        # mean free path and the definition still describe the study; the
        # number itself is left unset rather than invented.
        logger.debug("Knudsen number omitted for a firing with no swept "
                     "source distance (reference_length: source_distance)")
        return metadata
    metadata["knudsen_number"] = spec.knudsen_number(source_distance)
    metadata["knudsen_reference_length"] = spec.reference_length_for(
        source_distance)
    return metadata


def export_surface_distributions(config: StudyConfig,
                                 geometry: TargetGeometry, firing: Firing,
                                 per_face: dict[str, Any], *, case_id: str,
                                 firing_id: int, plate_angle_deg: float,
                                 source_distance: float,
                                 output_dir: str) -> dict[str, str]:
    """Write one panel-local distribution CSV per component, if enabled.

    An ADDITION to the VTK export, never a replacement: the same native
    per-face values, in the target's own (u, v) coordinates, with no
    interpolation. Returns ``{component name: csv path}``, empty when
    ``output.surface_distribution.enabled`` is false.
    """
    if not config.output.write_surface_distribution:
        return {}

    u_hat, v_hat, n_hat = config.target.local_basis()
    reference = config.target.reference_point
    knudsen = knudsen_metadata(config, source_distance)
    ensure_dir(output_dir)

    paths: dict[str, str] = {}
    for component_name, face_indices in geometry.components:
        rows = distribution_rows(
            face_indices, geometry.centroids, geometry.areas,
            per_face["pressures"], per_face["shear_stress"],
            per_face["heat_flux_rate"], per_face.get("strikes"),
            reference_point=reference, u_hat=u_hat, v_hat=v_hat)
        metadata: dict[str, Any] = {
            "study_name": config.study_name,
            "case_id": case_id,
            "component": component_name,
            "firing_id": firing_id,
            "geometry_id": config.target.geometry_id,
            "coordinate_system": config.coordinate_system,
            "panel_basis": {
                "origin": [float(value) for value in reference],
                "u": [float(value) for value in u_hat],
                "v": [float(value) for value in v_hat],
                "n": [float(value) for value in n_hat],
            },
            "plate_angle_deg": float(plate_angle_deg),
            "source_distance": float(source_distance),
            "source_offset_u": float(firing.source_offset_u),
            "source_offset_v": float(firing.source_offset_v),
            "source_axis_mode": firing.source_axis_mode,
            "plume_source_position": [float(v) for v in firing.position],
            "plume_model": config.plume_model,
            "plume_model_parameters": dict(config.plume_model_parameters),
            "firing_duration_s": float(firing.duration_s),
        }
        metadata.update(knudsen)

        name = f"{case_id}_{component_name}_firing{firing_id:03d}.csv"
        paths[component_name] = write_surface_distribution(
            os.path.join(output_dir, name), rows, metadata)
    return paths


def build_case_results(config: StudyConfig, geometry: TargetGeometry,
                       firing: Firing, per_face: dict[str, Any], *,
                       case_id: str, firing_id: int,
                       plate_angle_deg: float, source_distance: float,
                       jfh_path: str, vtk_path: str | None,
                       code_version_id: str, timestamp: str,
                       surface_distribution_paths: dict[str, str] | None = None,
                       ) -> list[CaseResult]:
    """Integrate one firing's per-face fields into one record per component.

    Each record also carries the resultants projected on the target's
    surface-local basis (normal force, local moments, panel-local center of
    pressure) and, when configured, the derived Knudsen metadata.
    """
    flow = flow_directions(geometry.centroids, firing.position)
    u_hat, v_hat, n_hat = config.target.local_basis()
    knudsen = knudsen_metadata(config, source_distance)
    distributions = surface_distribution_paths or {}
    records: list[CaseResult] = []
    for component_name, face_indices in geometry.components:
        loads = integrate_component_loads(
            component_name=component_name,
            face_indices=face_indices,
            centroids=geometry.centroids, unit_normals=geometry.normals,
            areas=geometry.areas,
            pressures=per_face["pressures"],
            shear_stresses=per_face["shear_stress"],
            heat_fluxes=per_face["heat_flux_rate"],
            strikes=per_face["strikes"],
            moment_reference_point=config.loads.moment_reference_point,
            flow_unit_vectors=flow,
            normalization=config.loads.normalization)

        records.append(CaseResult.from_loads(
            loads,
            panel=project_to_panel_frame(loads, u_hat, v_hat, n_hat),
            study_name=config.study_name,
            case_id=case_id,
            firing_id=firing_id,
            geometry_id=config.target.geometry_id,
            mesh_faces=geometry.n_faces,
            coordinate_system=config.coordinate_system,
            units=dict(config.units),
            plume_source_position=[float(v) for v in firing.position],
            plume_source_orientation=[
                float(v) for v in np.asarray(firing.dcm).ravel()],
            target_normal=[float(v) for v in config.target.normal],
            target_tangent=[float(v) for v in config.target.tangent],
            target_reference_point=[
                float(v) for v in config.target.reference_point],
            plate_angle_deg=float(plate_angle_deg),
            source_distance=float(source_distance),
            source_offset_u=float(firing.source_offset_u),
            source_offset_v=float(firing.source_offset_v),
            source_axis_mode=firing.source_axis_mode,
            firing_duration_s=float(firing.duration_s),
            thrusters=[int(t) for t in firing.thrusters],
            plume_model=config.plume_model,
            plume_model_parameters=dict(config.plume_model_parameters),
            vtk_path=vtk_path,
            jfh_path=jfh_path,
            surface_distribution_path=distributions.get(component_name),
            config_path=config.source_path,
            case_dir=os.path.abspath(config.case_dir),
            code_version=code_version_id,
            generated_at=timestamp,
            **knudsen))
    return records


def study_plots_for(config: StudyConfig, results: Any, out_dir: str,
                    comparison: Any = None) -> list[str]:
    """Generate whichever optional figures this study's sweep calls for.

    Both engines plot identically, so the choice lives here rather than in
    either of them (and never in the :class:`TradeStudy` façade):

    * the angle-sweep trends of :mod:`pyrpod.mdao.study_plots` are always
      produced, plus the study-vs-reference figure when a comparison exists;
    * an OFFSET sweep additionally gets the panel-local trends of
      :mod:`pyrpod.mdao.panel_plots`;
    * per-case panel pressure maps are drawn when
      ``output.plots.per_case_distribution`` is set AND the distributions
      they read were exported.

    Both plotting modules are imported lazily, so a study that asks for no
    plots never imports matplotlib.
    """
    from pyrpod.mdao import study_plots

    written = list(study_plots.plot_sweep_trends(results, out_dir,
                                                 comparison=comparison))

    sweeps_offsets = (config.sweep.source_axis_mode == "parallel_to_normal"
                      or len(config.sweep.source_offsets_u) > 1
                      or len(config.sweep.source_offsets_v) > 1)
    if not sweeps_offsets:
        return written

    from pyrpod.mdao import panel_plots

    written.extend(panel_plots.plot_offset_sweep_trends(results, out_dir))

    if config.output.write_distribution_plots:
        if not config.output.write_surface_distribution:
            logger.warning(
                "output.plots.per_case_distribution is set but "
                "output.surface_distribution.enabled is not; the per-case "
                "pressure maps are drawn from the exported distributions, "
                "so none were produced")
            return written
        seen: set[str] = set()
        for case in results.cases:
            if case.case_id in seen:
                continue
            seen.add(case.case_id)
            path = panel_plots.plot_panel_pressure_for_case(case, out_dir)
            if path:
                written.append(path)
    return written


def component_envelope(geometry: TargetGeometry,
                       cumulative: dict[str, Any],
                       ) -> dict[str, dict[str, float]]:
    """Worst-case per-face statistics over a whole firing history.

    Only meaningful when ONE history spans the sweep: the strike pipeline's
    cumulative arrays then hold the maximum pressure and shear any pose
    produced on each face, the summed heat-flux load, and the number of times
    each face was struck. Reported per component so a multi-component target
    keeps them separate.
    """
    envelope: dict[str, dict[str, float]] = {}
    for component_name, face_indices in geometry.components:
        indices = np.asarray(face_indices, dtype=np.int64)
        areas = geometry.areas[indices]
        strikes = np.asarray(cumulative["cum_strikes"], dtype=float)[indices]
        pressures = np.asarray(cumulative["max_pressures"],
                               dtype=float)[indices]
        shears = np.asarray(cumulative["max_shears"], dtype=float)[indices]
        heat_load = np.asarray(cumulative["cum_heat_flux_load"],
                               dtype=float)[indices]
        struck = strikes > 0.0
        envelope[component_name] = {
            "max_pressure": float(np.max(pressures)),
            "max_shear_stress": float(np.max(shears)),
            "max_heat_flux_load": float(np.max(heat_load)),
            "total_strike_events": float(np.sum(strikes)),
            "unique_struck_faces": int(np.count_nonzero(struck)),
            "swept_affected_area": float(np.sum(areas[struck])),
            "component_area": float(np.sum(areas)),
        }
    return envelope
