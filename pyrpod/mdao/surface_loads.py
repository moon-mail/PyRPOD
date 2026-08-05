"""
Integration of per-face plume-impingement fields into component-level loads.

The strike pipeline (:mod:`pyrpod.plume.PlumeStrikeCalculator`) produces
per-face pressure, shear-stress and heat-flux arrays. This module turns those
arrays into the integrated quantities an engineer reports for a target
component -- force, moment, center of pressure, peak surface loads, affected
area and (only when normalization inputs are supplied) nondimensional
coefficients. Per-face data are never discarded here; they keep flowing to
the VTK writer unchanged.

Sign and direction conventions
------------------------------
These follow the strike pipeline's own face-selection convention, in which a
face is struck only when its mesh unit normal ``n_hat`` opposes the plume
direction (``n_hat . u_plume < 0``), i.e. ``n_hat`` points back toward the
plume source:

* **Pressure** acts INTO the surface, along ``-n_hat``:
  ``F_p,i = -p_i * A_i * n_hat_i``.
* **Shear** acts along the tangential projection of the local flow direction
  ``u_hat_i`` (the collisionless flow is radial, so ``u_hat_i`` is the unit
  vector from the plume source to the face centroid):
  ``t_hat_i = normalize(u_hat_i - (u_hat_i . n_in_i) * n_in_i)``,
  ``n_in_i = -n_hat_i``, and ``F_s,i = tau_i * A_i * t_hat_i``.
  Faces at normal incidence have no tangential direction and contribute no
  shear force (their shear magnitude is zero there anyway).
* **Moment** about the user-supplied reference point ``r_ref``:
  ``M_ref = sum (r_i - r_ref) x F_i``.

Center of pressure
------------------
A unique three-dimensional center of pressure does not exist in general: the
resultant of a distributed load is a force plus a couple, and only the
component of the moment perpendicular to the force can be represented by
shifting the force's line of action. The definition used here is stated
explicitly:

    ``r_cop = r_ref + (F x M_ref) / |F|^2``

which is the point ON THE LINE OF ACTION of the resultant force that lies
CLOSEST to the moment reference point. The moment component parallel to the
force is reported separately as ``residual_couple`` -- it cannot be removed
by any choice of ``r_cop``. When the resultant force is zero or numerically
insignificant against the summed face-force magnitudes (near-total
cancellation, e.g. a closed target loaded symmetrically), no center of
pressure is returned and the status field records why; no misleadingly large
value is ever produced.

The auxiliary ``pressure_weighted_centroid``,
``sum(p_i A_i r_i) / sum(p_i A_i)``, is always defined for a non-zero
pressure load and coincides with the classical center of pressure for a
planar component under unidirectional pressure.

Panel-local reporting
---------------------
The global-frame resultants above are also reported on the target's own
surface-local basis ``(u_hat, v_hat, n_hat)``
(:meth:`pyrpod.mdao.study_config.TargetSpec.local_basis`), which is what a
panel study reads: a longitudinal axis ``u``, a transverse axis ``v``, and
the normal ``n``, with ``u x v = n``. :func:`project_to_panel_frame` builds
that projection, with these SIGN CONVENTIONS:

* ``normal_force = -F . n_hat``. The target normal points TOWARD the plume
  source, and the plume pushes into the surface, so a compressive
  impingement load is POSITIVE. A positive normal force therefore always
  means "pressed away from the source"; a negative one would mean the
  resultant pulls the panel toward the source.
* ``local_moment_u/v/n = M_ref . u_hat / v_hat / n_hat``, the components of
  the SAME moment vector the global-frame fields report, about the same
  reference point. Positive follows the right-hand rule about that axis.
  Worked sign: a pressure patch centred at ``+u`` pushes along ``-n``, so
  ``M = (u_off * u_hat) x (-F_n * n_hat) = +u_off * F_n * v_hat`` (using
  ``u x n = -v`` and ``F_n < 0`` along n). A source displaced toward ``+u``
  therefore produces a POSITIVE ``local_moment_v``, growing with the offset
  until the patch starts leaving the panel.
* ``center_of_pressure_u/v = (r_cop - r_ref) . u_hat / v_hat``, the
  panel-local coordinates of the center of pressure MEASURED FROM the
  moment reference point (the panel center for the ISS-panel studies).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from pyrpod.mdao.study_config import ComponentSpec, Normalization

__all__ = [
    "ComponentLoads",
    "PanelProjection",
    "face_areas",
    "flow_directions",
    "integrate_component_loads",
    "panel_local_coordinates",
    "project_to_panel_frame",
    "select_component_faces",
]

#: Relative floor below which a resultant force is treated as cancellation
#: noise rather than a meaningful resultant (see the module docstring).
FORCE_SIGNIFICANCE_RATIO = 1e-6


def face_areas(vectors: NDArray[np.float64]) -> NDArray[np.float64]:
    """Per-face triangle areas for an (N, 3, 3) array of face vertices."""
    vectors = np.asarray(vectors, dtype=float)
    v0, v1, v2 = vectors[:, 0], vectors[:, 1], vectors[:, 2]
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)


def flow_directions(centroids: NDArray[np.float64],
                    source_position: Sequence[float] | NDArray[np.float64],
                    ) -> NDArray[np.float64]:
    """Unit vectors from the plume source to each face centroid.

    The collisionless plume flow is radial from the nozzle exit, so this is
    the local flow direction at each face -- the same quantity the strike
    calculator uses for the true incidence angle.
    """
    centroids = np.asarray(centroids, dtype=float)
    relative = centroids - np.asarray(source_position, dtype=float).reshape(3)
    norms = np.linalg.norm(relative, axis=1)
    directions = np.zeros_like(relative)
    np.divide(relative, norms[:, None], out=directions,
              where=(norms > 0.0)[:, None])
    return directions


def select_component_faces(component: ComponentSpec,
                           centroids: NDArray[np.float64],
                           n_faces: int) -> NDArray[np.int64]:
    """Face indices belonging to one target component.

    Raises
    ------
    ValueError
        If an explicit face index is out of range, or a bounds selector
        matches no face (a silently empty component would produce
        meaningless zero loads).
    """
    if component.selector == "all":
        return np.arange(n_faces, dtype=np.int64)

    if component.selector == "face_indices":
        indices = np.asarray(component.face_indices, dtype=np.int64)
        if indices.min() < 0 or indices.max() >= n_faces:
            raise ValueError(
                f"component {component.name!r}: face index out of range for a "
                f"{n_faces}-face mesh")
        return indices

    lo = np.asarray(component.bounds["min"], dtype=float)
    hi = np.asarray(component.bounds["max"], dtype=float)
    inside = np.all((centroids >= lo) & (centroids <= hi), axis=1)
    indices = np.nonzero(inside)[0].astype(np.int64)
    if indices.size == 0:
        raise ValueError(
            f"component {component.name!r}: bounds selector matched no face "
            "of the target mesh")
    return indices


@dataclass(frozen=True)
class ComponentLoads:
    """Integrated loads for one target component and one firing.

    All vectors are in the case's global frame and SI units (N, N*m, m).
    ``coefficients`` is empty when the configuration supplies no (or
    incomplete) normalization inputs -- coefficients are never invented.
    """

    component: str
    n_faces: int
    n_struck_faces: int
    total_area: float
    affected_area: float

    pressure_force: NDArray[np.float64]
    shear_force: NDArray[np.float64]
    force: NDArray[np.float64]
    force_magnitude: float

    moment_reference_point: NDArray[np.float64]
    pressure_moment: NDArray[np.float64]
    shear_moment: NDArray[np.float64]
    moment: NDArray[np.float64]
    moment_magnitude: float

    center_of_pressure: NDArray[np.float64] | None
    center_of_pressure_status: str
    residual_couple: float
    pressure_weighted_centroid: NDArray[np.float64] | None

    max_pressure: float
    max_shear_stress: float
    max_heat_flux: float
    total_heat_load: float

    coefficients: dict[str, float] = field(default_factory=dict)

    @property
    def has_coefficients(self) -> bool:
        return bool(self.coefficients)


def integrate_component_loads(
    *,
    component_name: str,
    face_indices: NDArray[np.int64] | Sequence[int],
    centroids: NDArray[np.float64],
    unit_normals: NDArray[np.float64],
    areas: NDArray[np.float64],
    pressures: NDArray[np.float64],
    shear_stresses: NDArray[np.float64],
    heat_fluxes: NDArray[np.float64],
    strikes: NDArray[np.float64] | None,
    moment_reference_point: Sequence[float] | NDArray[np.float64],
    source_position: Sequence[float] | NDArray[np.float64] | None = None,
    flow_unit_vectors: NDArray[np.float64] | None = None,
    normalization: Normalization | None = None,
) -> ComponentLoads:
    """Integrate per-face fields into one component's resultant loads.

    Parameters
    ----------
    component_name : str
        Reported component identifier.
    face_indices : array-like of int
        Indices of the faces belonging to this component (see
        :func:`select_component_faces`).
    centroids, unit_normals, areas : np.ndarray
        Per-face geometry of the FULL target mesh, indexed by ``face_indices``.
    pressures, shear_stresses, heat_fluxes : np.ndarray
        Per-face load fields of the full mesh for one firing (Pa, Pa, W/m^2).
    strikes : np.ndarray or None
        Per-face strike counts, used for the affected-area tally. When None,
        faces carrying any non-zero load are counted instead.
    moment_reference_point : array-like
        User-defined point the moment is taken about.
    source_position : array-like, optional
        Plume-source position, used to derive the local flow directions the
        shear force acts along. Required unless ``flow_unit_vectors`` is given.
    flow_unit_vectors : np.ndarray, optional
        Precomputed per-face flow directions (see :func:`flow_directions`).
    normalization : Normalization, optional
        Coefficient normalization inputs. Coefficients are computed only for
        the ones whose inputs are complete.

    Returns
    -------
    ComponentLoads
    """
    indices = np.asarray(face_indices, dtype=np.int64)
    if indices.size == 0:
        raise ValueError(
            f"component {component_name!r} selects no faces; refusing to "
            "report zero loads for an empty component")

    centroid = np.asarray(centroids, dtype=float)[indices]
    normal = np.asarray(unit_normals, dtype=float)[indices]
    area = np.asarray(areas, dtype=float)[indices]
    pressure = np.asarray(pressures, dtype=float)[indices]
    shear = np.asarray(shear_stresses, dtype=float)[indices]
    heat_flux = np.asarray(heat_fluxes, dtype=float)[indices]

    if flow_unit_vectors is None:
        if source_position is None:
            raise ValueError(
                "integrate_component_loads needs either source_position or "
                "flow_unit_vectors to orient the shear force")
        flow = flow_directions(centroid, source_position)
    else:
        flow = np.asarray(flow_unit_vectors, dtype=float)[indices]

    # --- forces -----------------------------------------------------------
    # Pressure pushes into the surface (along -n_hat, the plume direction for
    # every struck face); shear drags along the tangential flow projection.
    pressure_face_force = -(pressure * area)[:, None] * normal

    inward = -normal
    cos_incidence = np.einsum("ij,ij->i", flow, inward)
    tangential = flow - cos_incidence[:, None] * inward
    tangential_norm = np.linalg.norm(tangential, axis=1)
    tangent_hat = np.zeros_like(tangential)
    np.divide(tangential, tangential_norm[:, None], out=tangent_hat,
              where=(tangential_norm > 1e-12)[:, None])
    shear_face_force = (shear * area)[:, None] * tangent_hat

    pressure_force = pressure_face_force.sum(axis=0)
    shear_force = shear_face_force.sum(axis=0)
    total_force = pressure_force + shear_force

    # --- moments ----------------------------------------------------------
    reference = np.asarray(moment_reference_point, dtype=float).reshape(3)
    arm = centroid - reference
    pressure_moment = np.cross(arm, pressure_face_force).sum(axis=0)
    shear_moment = np.cross(arm, shear_face_force).sum(axis=0)
    total_moment = pressure_moment + shear_moment

    # --- center of pressure ----------------------------------------------
    face_force_scale = float(np.sum(np.linalg.norm(
        pressure_face_force + shear_face_force, axis=1)))
    force_magnitude = float(np.linalg.norm(total_force))
    cop, cop_status, residual_couple = _center_of_pressure(
        total_force, total_moment, reference, force_magnitude,
        face_force_scale)

    pressure_weight = pressure * area
    weight_sum = float(np.sum(pressure_weight))
    pressure_centroid = (
        (pressure_weight[:, None] * centroid).sum(axis=0) / weight_sum
        if weight_sum > 0.0 else None)

    # --- surface-field statistics ----------------------------------------
    if strikes is None:
        affected = (pressure != 0.0) | (shear != 0.0) | (heat_flux != 0.0)
    else:
        affected = np.asarray(strikes, dtype=float)[indices] > 0.0

    loads = ComponentLoads(
        component=component_name,
        n_faces=int(indices.size),
        n_struck_faces=int(np.count_nonzero(affected)),
        total_area=float(np.sum(area)),
        affected_area=float(np.sum(area[affected])),
        pressure_force=pressure_force,
        shear_force=shear_force,
        force=total_force,
        force_magnitude=force_magnitude,
        moment_reference_point=reference,
        pressure_moment=pressure_moment,
        shear_moment=shear_moment,
        moment=total_moment,
        moment_magnitude=float(np.linalg.norm(total_moment)),
        center_of_pressure=cop,
        center_of_pressure_status=cop_status,
        residual_couple=residual_couple,
        pressure_weighted_centroid=pressure_centroid,
        max_pressure=float(np.max(pressure)) if pressure.size else 0.0,
        max_shear_stress=float(np.max(shear)) if shear.size else 0.0,
        max_heat_flux=float(np.max(heat_flux)) if heat_flux.size else 0.0,
        total_heat_load=float(np.sum(heat_flux * area)),
        coefficients={},
    )
    return _with_coefficients(loads, normalization)


@dataclass(frozen=True)
class PanelProjection:
    """Integrated loads projected onto a target's surface-local basis.

    See the module docstring for the sign conventions. Every field is a
    plain float (or None when the underlying quantity does not exist), so a
    result record carries them as flat, directly plottable columns.
    """

    normal_force: float
    local_force_u: float
    local_force_v: float
    local_moment_u: float
    local_moment_v: float
    local_moment_n: float
    center_of_pressure_u: float | None
    center_of_pressure_v: float | None


def panel_local_coordinates(
    points: NDArray[np.float64],
    reference_point: Sequence[float] | NDArray[np.float64],
    u_hat: Sequence[float] | NDArray[np.float64],
    v_hat: Sequence[float] | NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Panel-local ``(u, v)`` coordinates of points, measured from a reference.

    A pure projection of ``points - reference_point`` onto the two in-surface
    axes; nothing is interpolated, fitted or flattened. For a planar target
    whose points lie in the panel, ``(u, v)`` recovers the panel coordinates
    exactly.

    Parameters
    ----------
    points : np.ndarray
        (N, 3) points in the case's global frame (face centroids, typically).
    reference_point : array-like
        Panel-local origin, normally the target reference point.
    u_hat, v_hat : array-like
        The in-surface axes (see
        :meth:`pyrpod.mdao.study_config.TargetSpec.local_basis`).

    Returns
    -------
    (np.ndarray, np.ndarray)
        The ``u`` and ``v`` coordinate arrays, each of length N.
    """
    relative = (np.asarray(points, dtype=float)
                - np.asarray(reference_point, dtype=float).reshape(3))
    return (relative @ np.asarray(u_hat, dtype=float).reshape(3),
            relative @ np.asarray(v_hat, dtype=float).reshape(3))


def project_to_panel_frame(
    loads: ComponentLoads,
    u_hat: Sequence[float] | NDArray[np.float64],
    v_hat: Sequence[float] | NDArray[np.float64],
    n_hat: Sequence[float] | NDArray[np.float64],
) -> PanelProjection:
    """Project one component's resultants onto the surface-local basis.

    The projection is exact and adds no physics: it re-expresses the force
    and moment vectors this module already integrated, in the basis a panel
    study reports. See the module docstring for the sign conventions, in
    particular that ``normal_force`` is positive for a load pressing INTO
    the panel (away from the plume source).
    """
    u = np.asarray(u_hat, dtype=float).reshape(3)
    v = np.asarray(v_hat, dtype=float).reshape(3)
    n = np.asarray(n_hat, dtype=float).reshape(3)

    force = np.asarray(loads.force, dtype=float).reshape(3)
    moment = np.asarray(loads.moment, dtype=float).reshape(3)

    cop_u: float | None = None
    cop_v: float | None = None
    if loads.center_of_pressure is not None:
        arm = (np.asarray(loads.center_of_pressure, dtype=float).reshape(3)
               - np.asarray(loads.moment_reference_point,
                            dtype=float).reshape(3))
        cop_u = float(arm @ u)
        cop_v = float(arm @ v)

    return PanelProjection(
        # Positive into the panel: the normal points at the plume source.
        normal_force=float(-(force @ n)),
        local_force_u=float(force @ u),
        local_force_v=float(force @ v),
        local_moment_u=float(moment @ u),
        local_moment_v=float(moment @ v),
        local_moment_n=float(moment @ n),
        center_of_pressure_u=cop_u,
        center_of_pressure_v=cop_v,
    )


def _center_of_pressure(
    force: NDArray[np.float64], moment: NDArray[np.float64],
    reference: NDArray[np.float64], force_magnitude: float,
    face_force_scale: float,
) -> tuple[NDArray[np.float64] | None, str, float]:
    """Line-of-action center of pressure, or an explicit reason there is none.

    See the module docstring for the definition. Returns
    ``(point_or_None, status, residual_couple)`` where ``status`` is one of
    ``'ok'``, ``'zero_load'`` (nothing is loaded) or ``'ill_conditioned'``
    (the resultant is cancellation noise against the summed face forces).
    """
    if force_magnitude == 0.0:
        return None, "zero_load" if face_force_scale == 0.0 else \
            "ill_conditioned", float(np.linalg.norm(moment))
    if face_force_scale > 0.0 and (
            force_magnitude < FORCE_SIGNIFICANCE_RATIO * face_force_scale):
        return None, "ill_conditioned", float(np.linalg.norm(moment))

    direction = force / force_magnitude
    residual_couple = float(np.dot(moment, direction))
    point = reference + np.cross(force, moment) / force_magnitude ** 2
    return point, "ok", abs(residual_couple)


def _with_coefficients(loads: ComponentLoads,
                       normalization: Normalization | None) -> ComponentLoads:
    """Attach nondimensional coefficients, but only the fully normalizable ones."""
    if normalization is None:
        return loads

    coefficients: dict[str, float] = {}
    q_dyn = normalization.dynamic_pressure
    area_ref = normalization.reference_area
    length_ref = normalization.reference_length

    if q_dyn is not None:
        coefficients["Cp_max"] = loads.max_pressure / q_dyn
        coefficients["Cf_max"] = loads.max_shear_stress / q_dyn
    if normalization.reference_heat_flux is not None:
        coefficients["Cq_max"] = (loads.max_heat_flux
                                  / normalization.reference_heat_flux)
    if normalization.has_force_inputs:
        # mypy: has_force_inputs guarantees both are set
        denominator = float(q_dyn) * float(area_ref)  # type: ignore[arg-type]
        coefficients["CF"] = loads.force_magnitude / denominator
        coefficients["CFx"] = float(loads.force[0]) / denominator
        coefficients["CFy"] = float(loads.force[1]) / denominator
        coefficients["CFz"] = float(loads.force[2]) / denominator
    if normalization.has_moment_inputs:
        denominator = (float(q_dyn) * float(area_ref)  # type: ignore[arg-type]
                       * float(length_ref))            # type: ignore[arg-type]
        coefficients["CM"] = loads.moment_magnitude / denominator
        coefficients["CMx"] = float(loads.moment[0]) / denominator
        coefficients["CMy"] = float(loads.moment[1]) / denominator
        coefficients["CMz"] = float(loads.moment[2]) / denominator

    if not coefficients:
        return loads
    return ComponentLoads(**{**_as_kwargs(loads), "coefficients": coefficients})


def _as_kwargs(loads: ComponentLoads) -> Mapping[str, Any]:
    return {f: getattr(loads, f) for f in loads.__dataclass_fields__}
