"""
Structured results for prescribed plume/target validation studies.

One :class:`CaseResult` records everything needed to REPRODUCE and later
COMPARE a single (case, component, firing) evaluation: the study/case
identity, the geometry and mesh it ran on, the coordinate system and units,
the prescribed plume-source pose, the swept parameters, the plume model and
its parameters, the integrated loads, the surface-field peaks, the optional
coefficients, the VTK artifact path, and the code/configuration provenance.

Two machine-readable artifacts are written, both in formats already used by
the repository and with no new dependency:

* a flat CSV of every case row (one row per case x component x firing), for
  spreadsheets, plotting and quick diffing;
* a JSON metadata document carrying the study-level provenance plus the
  nested per-case records, which is the exchange format an externally
  generated dataset (DSMC, experiment, another code) is later transformed
  into for :mod:`pyrpod.mdao.reference_data`.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.typing import NDArray

from pyrpod.mdao.surface_loads import ComponentLoads, PanelProjection

__all__ = ["CaseResult", "StudyResults", "code_version"]


def code_version(repo_dir: str | os.PathLike[str] | None = None) -> str:
    """Short git commit of the working tree, or ``'unknown'``.

    Provenance is best-effort: a study must still run in an exported
    source tree with no git metadata available.
    """
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.fspath(repo_dir) if repo_dir else None,
            capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    if completed.returncode != 0:
        return "unknown"
    return completed.stdout.strip() or "unknown"


def _vector(value: Sequence[float] | NDArray[np.float64] | None
            ) -> list[float] | None:
    if value is None:
        return None
    return [float(v) for v in np.asarray(value, dtype=float).reshape(-1)]


@dataclass
class CaseResult:
    """Result record for one case, component and firing.

    Attributes are grouped as: identity, geometry, pose/sweep parameters,
    model, integrated loads, surface-field statistics, coefficients, and
    artifacts/provenance. Every field is plain data (str, float, list) so the
    record serializes to CSV and JSON without custom encoders.
    """

    # --- identity ---------------------------------------------------------
    study_name: str
    case_id: str
    component: str
    firing_id: int

    # --- geometry ---------------------------------------------------------
    geometry_id: str
    mesh_faces: int
    component_faces: int
    component_area: float
    coordinate_system: str
    units: dict[str, str]

    # --- prescribed pose and swept parameters -----------------------------
    plume_source_position: list[float]
    plume_source_orientation: list[float]
    target_normal: list[float]
    target_tangent: list[float]
    target_reference_point: list[float]
    plate_angle_deg: float
    source_distance: float
    firing_duration_s: float
    thrusters: list[int]

    # --- model ------------------------------------------------------------
    plume_model: str
    plume_model_parameters: dict[str, Any]

    # --- integrated loads -------------------------------------------------
    pressure_force: list[float]
    shear_force: list[float]
    force: list[float]
    force_magnitude: float
    moment_reference_point: list[float]
    pressure_moment: list[float]
    shear_moment: list[float]
    moment: list[float]
    moment_magnitude: float
    center_of_pressure: list[float] | None
    center_of_pressure_status: str
    residual_couple: float
    pressure_weighted_centroid: list[float] | None

    # --- surface-field statistics ----------------------------------------
    max_pressure: float
    max_shear_stress: float
    max_heat_flux: float
    total_heat_load: float
    affected_area: float
    struck_faces: int

    # --- optional coefficients -------------------------------------------
    coefficients: dict[str, float] = field(default_factory=dict)
    coefficients_available: bool = False

    # --- panel-local pose (offset sweeps) --------------------------------
    # Zero / 'aim_at_reference' for every study that does not sweep offsets,
    # so an older result file round-trips unchanged.
    source_offset_u: float = 0.0
    source_offset_v: float = 0.0
    source_axis_mode: str = "aim_at_reference"

    # --- panel-local integrated loads ------------------------------------
    # See pyrpod.mdao.surface_loads for the sign conventions: normal_force is
    # POSITIVE for a load pressing into the panel (away from the source);
    # the moments are the global moment vector projected on (u, v, n) about
    # the same reference point; the centers of pressure are measured FROM
    # that reference point along u and v.
    normal_force: float | None = None
    local_force_u: float | None = None
    local_force_v: float | None = None
    local_moment_u: float | None = None
    local_moment_v: float | None = None
    local_moment_n: float | None = None
    center_of_pressure_u: float | None = None
    center_of_pressure_v: float | None = None

    # --- derived Knudsen metadata (never an input to the physics) --------
    knudsen_number: float | None = None
    mean_free_path: float | None = None
    knudsen_reference_length: float | None = None
    knudsen_definition: str | None = None

    # --- artifacts and provenance ----------------------------------------
    vtk_path: str | None = None
    jfh_path: str | None = None
    surface_distribution_path: str | None = None
    config_path: str = ""
    case_dir: str = ""
    code_version: str = "unknown"
    generated_at: str = ""

    # ------------------------------------------------------------ accessors
    @property
    def model_variant(self) -> str:
        """Short label of the plume model that produced this record.

        ``'Simplified'`` or ``'Collisionless'`` -- the ``GasKinetics`` suffix
        of :attr:`plume_model` carries no information and makes a plot legend
        or a grouped CSV column needlessly wide.
        """
        return self.plume_model.replace("GasKinetics", "") or self.plume_model

    # ------------------------------------------------------------ builders
    @classmethod
    def from_loads(cls, loads: ComponentLoads,
                   panel: PanelProjection | None = None,
                   **metadata: Any) -> "CaseResult":
        """Build a result record from integrated loads plus study metadata.

        Parameters
        ----------
        loads : ComponentLoads
            Integrated global-frame loads for one component and firing.
        panel : PanelProjection, optional
            The same resultants projected on the target's surface-local
            basis (see :func:`pyrpod.mdao.surface_loads.project_to_panel_frame`).
            When omitted every panel-local field stays None, so a caller that
            has no surface basis is never given invented numbers.
        """
        panel_fields: dict[str, Any] = {}
        if panel is not None:
            panel_fields = {
                "normal_force": panel.normal_force,
                "local_force_u": panel.local_force_u,
                "local_force_v": panel.local_force_v,
                "local_moment_u": panel.local_moment_u,
                "local_moment_v": panel.local_moment_v,
                "local_moment_n": panel.local_moment_n,
                "center_of_pressure_u": panel.center_of_pressure_u,
                "center_of_pressure_v": panel.center_of_pressure_v,
            }
        return cls(
            **panel_fields,
            component=loads.component,
            component_faces=loads.n_faces,
            component_area=loads.total_area,
            pressure_force=_vector(loads.pressure_force) or [],
            shear_force=_vector(loads.shear_force) or [],
            force=_vector(loads.force) or [],
            force_magnitude=loads.force_magnitude,
            moment_reference_point=_vector(loads.moment_reference_point) or [],
            pressure_moment=_vector(loads.pressure_moment) or [],
            shear_moment=_vector(loads.shear_moment) or [],
            moment=_vector(loads.moment) or [],
            moment_magnitude=loads.moment_magnitude,
            center_of_pressure=_vector(loads.center_of_pressure),
            center_of_pressure_status=loads.center_of_pressure_status,
            residual_couple=loads.residual_couple,
            pressure_weighted_centroid=_vector(
                loads.pressure_weighted_centroid),
            max_pressure=loads.max_pressure,
            max_shear_stress=loads.max_shear_stress,
            max_heat_flux=loads.max_heat_flux,
            total_heat_load=loads.total_heat_load,
            affected_area=loads.affected_area,
            struck_faces=loads.n_struck_faces,
            coefficients=dict(loads.coefficients),
            coefficients_available=loads.has_coefficients,
            **metadata,
        )

    # --------------------------------------------------------- serialization
    def to_dict(self) -> dict[str, Any]:
        """Nested dictionary form (JSON metadata document)."""
        return asdict(self)

    def to_row(self) -> dict[str, Any]:
        """Flat dictionary form (CSV row).

        Vectors are expanded to ``<name>_x/_y/_z`` columns and coefficients to
        one column each, so the CSV is directly plottable.
        """
        row: dict[str, Any] = {
            "study_name": self.study_name,
            "case_id": self.case_id,
            "component": self.component,
            "firing_id": self.firing_id,
            "geometry_id": self.geometry_id,
            "mesh_faces": self.mesh_faces,
            "component_faces": self.component_faces,
            "component_area": self.component_area,
            "coordinate_system": self.coordinate_system,
            "plate_angle_deg": self.plate_angle_deg,
            "source_distance": self.source_distance,
            "source_offset_u": self.source_offset_u,
            "source_offset_v": self.source_offset_v,
            "source_axis_mode": self.source_axis_mode,
            "firing_duration_s": self.firing_duration_s,
            "thrusters": " ".join(str(t) for t in self.thrusters),
            "plume_model": self.plume_model,
            "model_variant": self.model_variant,
        }
        _expand(row, "plume_source_position", self.plume_source_position)
        _expand(row, "target_normal", self.target_normal)
        _expand(row, "target_tangent", self.target_tangent)
        _expand(row, "moment_reference_point", self.moment_reference_point)
        _expand(row, "pressure_force", self.pressure_force)
        _expand(row, "shear_force", self.shear_force)
        _expand(row, "force", self.force)
        row["force_magnitude"] = self.force_magnitude
        _expand(row, "pressure_moment", self.pressure_moment)
        _expand(row, "shear_moment", self.shear_moment)
        _expand(row, "moment", self.moment)
        row["moment_magnitude"] = self.moment_magnitude
        _expand(row, "center_of_pressure", self.center_of_pressure)
        row["center_of_pressure_status"] = self.center_of_pressure_status
        row["residual_couple"] = self.residual_couple
        _expand(row, "pressure_weighted_centroid",
                self.pressure_weighted_centroid)
        # Panel-local resultants: flat, directly plottable columns. Empty
        # rather than zero when the study supplied no surface basis.
        for name in ("normal_force", "local_force_u", "local_force_v",
                     "local_moment_u", "local_moment_v", "local_moment_n",
                     "center_of_pressure_u", "center_of_pressure_v"):
            value = getattr(self, name)
            row[name] = "" if value is None else float(value)
        row.update({
            "max_pressure": self.max_pressure,
            "max_shear_stress": self.max_shear_stress,
            "max_heat_flux": self.max_heat_flux,
            "total_heat_load": self.total_heat_load,
            "affected_area": self.affected_area,
            "struck_faces": self.struck_faces,
            "coefficients_available": self.coefficients_available,
        })
        for name, value in sorted(self.coefficients.items()):
            row[f"coeff_{name}"] = value
        # Derived Knudsen metadata; empty columns when the study configured
        # no knudsen block, never a fabricated value.
        for name in ("knudsen_number", "mean_free_path",
                     "knudsen_reference_length"):
            value = getattr(self, name)
            row[name] = "" if value is None else float(value)
        row["knudsen_definition"] = self.knudsen_definition or ""
        row.update({
            "vtk_path": self.vtk_path or "",
            "jfh_path": self.jfh_path or "",
            "surface_distribution_path": self.surface_distribution_path or "",
            "config_path": self.config_path,
            "case_dir": self.case_dir,
            "code_version": self.code_version,
            "generated_at": self.generated_at,
        })
        return row

    # ------------------------------------------------------------ accessors
    def quantity(self, name: str) -> float | list[float] | None:
        """Look up a comparable quantity by name (see reference_data).

        Supports the integrated loads, peaks, areas and coefficients under
        stable names, returning None when the quantity is unavailable for
        this result (e.g. a coefficient with no normalization inputs).
        """
        direct: dict[str, float | list[float] | None] = {
            "force": self.force,
            "force_magnitude": self.force_magnitude,
            "pressure_force": self.pressure_force,
            "shear_force": self.shear_force,
            "moment": self.moment,
            "moment_magnitude": self.moment_magnitude,
            "center_of_pressure": self.center_of_pressure,
            "max_pressure": self.max_pressure,
            "max_shear_stress": self.max_shear_stress,
            "max_heat_flux": self.max_heat_flux,
            "total_heat_load": self.total_heat_load,
            "affected_area": self.affected_area,
            # Panel-local scalars; None when no surface basis was supplied,
            # which quantity() already reports as "unavailable".
            "normal_force": self.normal_force,
            "local_force_u": self.local_force_u,
            "local_force_v": self.local_force_v,
            "local_moment_u": self.local_moment_u,
            "local_moment_v": self.local_moment_v,
            "local_moment_n": self.local_moment_n,
            "center_of_pressure_u": self.center_of_pressure_u,
            "center_of_pressure_v": self.center_of_pressure_v,
            # Derived metadata, comparable across cases like any other
            # scalar; it never participated in producing the loads.
            "knudsen_number": self.knudsen_number,
        }
        if name in direct:
            return direct[name]
        if name in self.coefficients:
            return self.coefficients[name]
        return None


def _expand(row: dict[str, Any], name: str,
            vector: Sequence[float] | None) -> None:
    for axis, index in (("x", 0), ("y", 1), ("z", 2)):
        row[f"{name}_{axis}"] = (float(vector[index])
                                 if vector is not None else "")


@dataclass
class StudyResults:
    """All case results of one study run, plus study-level provenance."""

    study_name: str
    output_dir: str
    provenance: dict[str, Any]
    cases: list[CaseResult] = field(default_factory=list)
    summary_csv_path: str | None = None
    metadata_path: str | None = None
    plot_paths: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.cases)

    def __iter__(self) -> Iterable[CaseResult]:  # type: ignore[override]
        return iter(self.cases)

    def rows(self) -> list[dict[str, Any]]:
        return [case.to_row() for case in self.cases]

    def for_component(self, component: str) -> list[CaseResult]:
        return [case for case in self.cases if case.component == component]

    # --------------------------------------------------------------- output
    def write_csv(self, path: str | os.PathLike[str]) -> str:
        """Write the flat per-case summary CSV; returns the path written."""
        rows = self.rows()
        if not rows:
            raise ValueError("no case results to write")
        columns: list[str] = []
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)
        path = os.fspath(path)
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        self.summary_csv_path = path
        return path

    def write_metadata(self, path: str | os.PathLike[str]) -> str:
        """Write the structured JSON metadata document; returns the path."""
        path = os.fspath(path)
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        document = {
            "schema": "pyrpod.plume_validation_study/1",
            "study_name": self.study_name,
            "generated_at": datetime.now(timezone.utc).isoformat(
                timespec="seconds"),
            "provenance": self.provenance,
            "n_cases": len(self.cases),
            "cases": [case.to_dict() for case in self.cases],
        }
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(document, handle, indent=2, sort_keys=False)
            handle.write("\n")
        self.metadata_path = path
        return path
