"""
YAML study configuration for prescribed plume/target validation sweeps.

A study configuration is a thin, declarative layer ON TOP OF an existing
PyRPOD case directory: the case's ``config.ini`` keeps owning the vehicle,
thruster, plume-model and target assets (so every existing case and public
API keeps working unchanged), while the YAML file adds only what a trade
study needs and a ``config.ini`` cannot express -- swept plate angles and
source distances, prescribed firing counts, the moment reference point,
coefficient normalization, and output/plot/reference settings.

The schema is deliberately explicit. Nothing is inferred that the caller
did not write down:

* the plume model is named explicitly and must be one of the collisionless
  models in :data:`pyrpod.plume.gas_kinetics_models.PLUME_MODELS`; the name
  SELECTS the model that computes the plume field, it is not just metadata;
* coefficients are computed only when every normalization input is present
  (see :class:`Normalization`); otherwise they are reported as unavailable;
* ``n_firings`` means the exact number of entries written to the Jet Firing
  History -- a mismatch against an explicitly prescribed firing list is an
  error, never a silent truncation;
* the Knudsen number is DERIVED METADATA only (see :class:`KnudsenSpec`):
  the mean free path must be supplied, is never inferred from gas
  properties, and never changes the analytical plume solution.

Example
-------
>>> cfg = StudyConfig.from_yaml('case/.../flat_plate_baseline.yaml')
>>> cfg.study_name
'cai2016_flat_plate_baseline'
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from typing import Any, Mapping, Sequence

import numpy as np
import yaml
from numpy.typing import NDArray

from pyrpod.plume.gas_kinetics_models import (
    DEFAULT_PLUME_MODEL,
    PLUME_MODELS,
    PlumeModelError,
    resolve_model_name,
)
from pyrpod.rpod.approach_maneuvers import (
    validate_n_firings as _validate_n_firings,
)

#: Plume model assumed when a configuration names none. Every configuration
#: written before model dispatch existed therefore keeps its exact behavior.
SUPPORTED_PLUME_MODEL = DEFAULT_PLUME_MODEL

#: Collisionless plume models a study may select, by class name.
SUPPORTED_PLUME_MODELS: tuple[str, ...] = tuple(sorted(PLUME_MODELS))

#: How the plume axis is oriented at each generated pose.
#:
#: ``aim_at_reference``
#:     The historical (and default) behavior: the source sits on an arc of
#:     radius ``source_distance`` about the target reference point and its
#:     axis points back at that point, so a swept ``plate_angles_deg``
#:     changes the approach angle.
#: ``parallel_to_normal``
#:     The source is TRANSLATED parallel to the target surface and its axis
#:     stays anti-parallel to the target normal, so the plume centerline
#:     strikes the surface at the requested panel-local offset. This is the
#:     ISS-panel convention; see :class:`SweepSpec`.
SOURCE_AXIS_MODES = ("aim_at_reference", "parallel_to_normal")

#: Default unit labels carried into the result metadata. PyRPOD works in SI
#: throughout; recording them makes an exported result self-describing.
DEFAULT_UNITS: dict[str, str] = {
    "length": "m",
    "area": "m^2",
    "force": "N",
    "moment": "N*m",
    "pressure": "Pa",
    "shear_stress": "Pa",
    "heat_flux": "W/m^2",
    "time": "s",
    "angle": "deg",
}


class StudyConfigError(ValueError):
    """Raised when a study configuration is missing or internally inconsistent."""


def _require(mapping: Mapping[str, Any], key: str, context: str) -> Any:
    if key not in mapping:
        raise StudyConfigError(
            f"missing required key {key!r} in {context}")
    return mapping[key]


def _as_vector(value: Any, context: str) -> NDArray[np.float64]:
    """Parse a 3-element vector, failing loudly on anything else."""
    try:
        vector = np.asarray(value, dtype=float).reshape(-1)
    except (TypeError, ValueError) as exc:
        raise StudyConfigError(
            f"{context}: expected three numbers, got {value!r}") from exc
    if vector.size != 3:
        raise StudyConfigError(
            f"{context}: expected three numbers, got {value!r}")
    return vector


def _as_float_list(value: Any, context: str) -> list[float]:
    if value is None or isinstance(value, (str, bytes)) or not isinstance(
            value, Sequence):
        raise StudyConfigError(f"{context}: expected a list of numbers, "
                               f"got {value!r}")
    try:
        return [float(v) for v in value]
    except (TypeError, ValueError) as exc:
        raise StudyConfigError(f"{context}: expected a list of numbers, "
                               f"got {value!r}") from exc


def _as_offsets(value: Any, context: str) -> list[float]:
    """Parse an optional panel-local offset list, defaulting to ``[0.0]``.

    An omitted (or explicitly null) list means "no offset", which is the
    behavior of every configuration written before offsets existed. Non-finite
    values are rejected rather than propagated into a pose.
    """
    if value is None:
        return [0.0]
    offsets = _as_float_list(value, context)
    if not offsets:
        raise StudyConfigError(
            f"{context} must not be empty; omit it for the default [0.0]")
    for offset in offsets:
        if not np.isfinite(offset):
            raise StudyConfigError(
                f"{context}: offsets must be finite, got {offset!r}")
    return offsets


@dataclass(frozen=True)
class ComponentSpec:
    """One target component: a named subset of the target mesh's faces.

    A component is selected either by ``face_indices`` (explicit) or by
    ``bounds`` (an axis-aligned box on face centroids, in the case's global
    frame). ``selector: all`` -- the default -- takes the whole mesh, which
    is the right answer for the single-plate and single-cylinder targets this
    branch ships. Nothing here assumes a flat plate.

    Attributes
    ----------
    name : str
        Component identifier, reported with every result row.
    selector : str
        ``'all'``, ``'face_indices'`` or ``'bounds'``.
    face_indices : tuple of int
        Explicit face indices when ``selector == 'face_indices'``.
    bounds : dict
        ``{'min': [x, y, z], 'max': [x, y, z]}`` when ``selector == 'bounds'``.
    """

    name: str
    selector: str = "all"
    face_indices: tuple[int, ...] = ()
    bounds: dict[str, tuple[float, float, float]] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ComponentSpec":
        name = str(_require(data, "name", "target.components entry"))
        if "face_indices" in data and data["face_indices"] is not None:
            indices = tuple(int(i) for i in data["face_indices"])
            if not indices:
                raise StudyConfigError(
                    f"target component {name!r}: face_indices is empty")
            return cls(name=name, selector="face_indices", face_indices=indices)
        if "bounds" in data and data["bounds"] is not None:
            raw = data["bounds"]
            lo = _as_vector(_require(raw, "min", f"component {name!r} bounds"),
                            f"component {name!r} bounds.min")
            hi = _as_vector(_require(raw, "max", f"component {name!r} bounds"),
                            f"component {name!r} bounds.max")
            if np.any(hi < lo):
                raise StudyConfigError(
                    f"target component {name!r}: bounds.max must be >= "
                    "bounds.min componentwise")
            bounds = {"min": (float(lo[0]), float(lo[1]), float(lo[2])),
                      "max": (float(hi[0]), float(hi[1]), float(hi[2]))}
            return cls(name=name, selector="bounds", bounds=bounds)
        selector = str(data.get("selector", "all"))
        if selector != "all":
            raise StudyConfigError(
                f"target component {name!r}: selector {selector!r} needs "
                "'face_indices' or 'bounds' to be supplied")
        return cls(name=name, selector="all")


@dataclass(frozen=True)
class TargetSpec:
    """Target geometry description and its component breakdown.

    ``reference_point`` is the geometric reference the sweep is built about
    (the plate center, the cylinder axis midpoint, ...). ``normal`` and
    ``tangent`` define the plane the approach angle is swept in: angle 0
    places the plume source on ``normal``, positive angles rotate it toward
    ``tangent``. They are geometry properties, not physics, so a curved
    target simply supplies the axes its sweep should use.

    The same two vectors also define the SURFACE-LOCAL basis used by the
    panel-local offsets, moments, center of pressure and distribution
    exports (see :meth:`local_basis`); no separate axis keys are needed.
    """

    geometry_id: str
    reference_point: NDArray[np.float64]
    normal: NDArray[np.float64]
    tangent: NDArray[np.float64]
    components: tuple[ComponentSpec, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any],
                     default_geometry_id: str) -> "TargetSpec":
        normal = _as_vector(data.get("normal", [0.0, 0.0, 1.0]),
                            "target.normal")
        tangent = _as_vector(data.get("tangent", [1.0, 0.0, 0.0]),
                             "target.tangent")
        for label, vec in (("normal", normal), ("tangent", tangent)):
            if not np.isfinite(vec).all() or np.linalg.norm(vec) == 0.0:
                raise StudyConfigError(
                    f"target.{label} must be a finite non-zero vector")
        normal = normal / np.linalg.norm(normal)
        tangent = tangent / np.linalg.norm(tangent)
        if abs(float(normal @ tangent)) > 1e-8:
            raise StudyConfigError(
                "target.normal and target.tangent must be orthogonal "
                f"(dot product {float(normal @ tangent):.3g})")

        raw_components = data.get("components")
        if raw_components is None:
            components: tuple[ComponentSpec, ...] = (
                ComponentSpec(name="target"),)
        else:
            if not isinstance(raw_components, Sequence) or not raw_components:
                raise StudyConfigError(
                    "target.components must be a non-empty list")
            components = tuple(ComponentSpec.from_mapping(entry)
                               for entry in raw_components)
            names = [component.name for component in components]
            if len(set(names)) != len(names):
                raise StudyConfigError(
                    f"target.components names must be unique, got {names}")

        return cls(
            geometry_id=str(data.get("geometry_id") or default_geometry_id),
            reference_point=_as_vector(data.get("reference_point",
                                                [0.0, 0.0, 0.0]),
                                       "target.reference_point"),
            normal=normal,
            tangent=tangent,
            components=components,
        )

    def local_basis(self) -> tuple[NDArray[np.float64], NDArray[np.float64],
                                   NDArray[np.float64]]:
        """Right-handed surface-local basis ``(u_hat, v_hat, n_hat)``.

        * ``u_hat`` is the target tangent -- the LONGITUDINAL in-surface
          axis, the 22 m dimension of the ISS-representative panel;
        * ``v_hat = n_hat x u_hat`` is the TRANSVERSE in-surface axis, the
          12 m dimension;
        * ``n_hat`` is the target normal, pointing toward the plume source.

        The triad satisfies ``u x v = n``, and ``v_hat`` is exactly the
        binormal :func:`pyrpod.mdao.firing_plan.pose_for` already uses for
        the second column of its DCM, so the pose convention and the
        panel-local reporting convention are the same basis.
        """
        n_hat = self.normal / np.linalg.norm(self.normal)
        u_hat = self.tangent / np.linalg.norm(self.tangent)
        v_hat = np.cross(n_hat, u_hat)
        return u_hat, v_hat / np.linalg.norm(v_hat), n_hat


@dataclass(frozen=True)
class PrescribedFiringSpec:
    """One explicitly prescribed firing: pose, active thrusters, duration."""

    position: NDArray[np.float64]
    dcm: NDArray[np.float64]
    thrusters: tuple[int, ...]
    duration_s: float

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], index: int,
                     default_thrusters: Sequence[int],
                     default_duration_s: float) -> "PrescribedFiringSpec":
        context = f"firings[{index}]"
        position = _as_vector(_require(data, "position", context),
                              f"{context}.position")
        dcm = np.asarray(_require(data, "dcm", context), dtype=float)
        if dcm.shape != (3, 3):
            raise StudyConfigError(
                f"{context}.dcm must be a 3x3 matrix, got shape {dcm.shape}")
        thrusters = tuple(int(t) for t in
                          data.get("thrusters", default_thrusters))
        if not thrusters:
            raise StudyConfigError(f"{context}.thrusters must not be empty")
        duration = float(data.get("duration_s", default_duration_s))
        if duration <= 0.0:
            raise StudyConfigError(
                f"{context}.duration_s must be positive, got {duration}")
        return cls(position=position, dcm=dcm, thrusters=thrusters,
                   duration_s=duration)


#: Sweep execution modes. ``per_case`` writes one Jet Firing History per
#: angle-distance combination (isolated poses); ``single_jfh`` writes ONE
#: history spanning the whole sweep, so every pose's strikes come from a
#: single pipeline run and the cumulative fields form a sweep envelope.
SWEEP_MODES = ("per_case", "single_jfh")


@dataclass(frozen=True)
class SweepPose:
    """One swept plume-source pose: angle, distance and panel-local offsets.

    The full parameterization of a generated pose, in the order
    :meth:`SweepSpec.sweep_poses` enumerates them. ``source_offset_u`` and
    ``source_offset_v`` are zero for every configuration written before the
    offset sweep existed, so such a study's poses are exactly what they
    always were.
    """

    plate_angle_deg: float
    source_distance: float
    source_offset_u: float = 0.0
    source_offset_v: float = 0.0
    axis_mode: str = "aim_at_reference"

    @property
    def key(self) -> tuple[float, float]:
        """The legacy ``(angle, distance)`` pose key."""
        return (self.plate_angle_deg, self.source_distance)


@dataclass(frozen=True)
class SweepSpec:
    """Parameter sweep and firing-count definition.

    Attributes
    ----------
    plate_angles_deg : tuple of float
        Approach angles swept in the target's (normal, tangent) plane;
        0 deg is head-on along the target normal. Meaningful only in
        ``aim_at_reference`` axis mode.
    source_distances : tuple of float
        Plume-source distances from ``TargetSpec.reference_point``, in the
        case's length units.
    source_offsets_u, source_offsets_v : tuple of float
        Panel-local translations of the plume source along the target's
        longitudinal (``u``) and transverse (``v``) axes (see
        :meth:`TargetSpec.local_basis`). Both default to ``(0.0,)``, which
        reproduces the on-axis poses exactly. Requires
        ``source_axis_mode: parallel_to_normal``.
    source_axis_mode : str
        ``'aim_at_reference'`` (default) or ``'parallel_to_normal'``; see
        :data:`SOURCE_AXIS_MODES`.
    n_firings : int
        Number of Jet Firing History entries contributed by EACH pose. In
        ``per_case`` mode that is the exact length of every case's history;
        in ``single_jfh`` mode the one shared history holds exactly
        ``len(poses) * n_firings`` entries. Either way the count is exact.
    mode : str
        ``'per_case'`` (default) or ``'single_jfh'``; see :data:`SWEEP_MODES`.
    firing_duration_s : float
        Firing time recorded for each JFH entry.
    thrusters : tuple of int
        JFH thruster indices active in every generated firing.
    firings : tuple of PrescribedFiringSpec
        Explicitly prescribed firings, replacing the generated poses. In
        ``per_case`` mode their count must equal ``n_firings``; in
        ``single_jfh`` mode it must equal ``len(poses) * n_firings``, and
        they are assigned to the poses in sweep order.
    """

    plate_angles_deg: tuple[float, ...]
    source_distances: tuple[float, ...]
    n_firings: int
    firing_duration_s: float
    thrusters: tuple[int, ...]
    firings: tuple[PrescribedFiringSpec, ...] = ()
    mode: str = "per_case"
    source_offsets_u: tuple[float, ...] = (0.0,)
    source_offsets_v: tuple[float, ...] = (0.0,)
    source_axis_mode: str = "aim_at_reference"

    @property
    def sweep_poses(self) -> tuple[SweepPose, ...]:
        """Fully parameterized poses in execution order.

        Distance-major, matching the committed sweep-JFH generators, then
        u offset, then v offset, then approach angle:

            for distance: for u_offset: for v_offset: for angle

        With the default single-element offset lists this collapses to the
        historical "all angles at the first distance, then all angles at the
        next", so an existing configuration's pose order is unchanged. For
        an offset sweep (a single angle, as ``parallel_to_normal`` requires)
        the length is exactly
        ``n_distances * n_u_offsets * n_v_offsets``.
        """
        return tuple(
            SweepPose(plate_angle_deg=angle, source_distance=distance,
                      source_offset_u=u_offset, source_offset_v=v_offset,
                      axis_mode=self.source_axis_mode)
            for distance in self.source_distances
            for u_offset in self.source_offsets_u
            for v_offset in self.source_offsets_v
            for angle in self.plate_angles_deg)

    @property
    def poses(self) -> tuple[tuple[float, float], ...]:
        """(plate angle, source distance) pairs in execution order.

        The legacy projection of :attr:`sweep_poses`, kept unchanged for
        callers that only key on the swept angle and distance. It has the
        same length and order as ``sweep_poses``; when offsets are swept,
        several entries share an ``(angle, distance)`` pair.
        """
        return tuple(pose.key for pose in self.sweep_poses)

    @property
    def total_firings(self) -> int:
        """Total JFH entries across the whole sweep."""
        return len(self.sweep_poses) * self.n_firings

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "SweepSpec":
        angles = tuple(_as_float_list(
            data.get("plate_angles_deg", [0.0]), "sweep.plate_angles_deg"))
        distances = tuple(_as_float_list(
            data.get("source_distances", []), "sweep.source_distances"))
        if not angles:
            raise StudyConfigError("sweep.plate_angles_deg must not be empty")
        if not distances:
            raise StudyConfigError("sweep.source_distances must not be empty")
        if any(d <= 0.0 for d in distances):
            raise StudyConfigError(
                "sweep.source_distances must all be positive")

        offsets_u = tuple(_as_offsets(data.get("source_offsets_u"),
                                      "sweep.source_offsets_u"))
        offsets_v = tuple(_as_offsets(data.get("source_offsets_v"),
                                      "sweep.source_offsets_v"))

        axis_mode = str(data.get("source_axis_mode", "aim_at_reference"))
        if axis_mode not in SOURCE_AXIS_MODES:
            raise StudyConfigError(
                f"sweep.source_axis_mode must be one of "
                f"{list(SOURCE_AXIS_MODES)}, got {axis_mode!r}")

        # Pose definitions that mean different things are never silently
        # combined: an aimed arc has no panel-local offset, and a fixed
        # axis parallel to -n has no approach angle.
        offsets_swept = (offsets_u != (0.0,) or offsets_v != (0.0,))
        if axis_mode == "aim_at_reference" and offsets_swept:
            raise StudyConfigError(
                "sweep.source_offsets_u / source_offsets_v translate the "
                "plume source parallel to the target surface, which is only "
                "defined for sweep.source_axis_mode: 'parallel_to_normal'; "
                f"got {axis_mode!r}")
        if axis_mode == "parallel_to_normal" and angles != (0.0,):
            raise StudyConfigError(
                "sweep.source_axis_mode: 'parallel_to_normal' fixes the "
                "plume axis anti-parallel to the target normal, so an "
                "approach angle has no meaning; remove "
                f"sweep.plate_angles_deg (got {list(angles)}) and sweep "
                "source_offsets_u / source_offsets_v instead")

        mode = str(data.get("mode", "per_case"))
        if mode not in SWEEP_MODES:
            raise StudyConfigError(
                f"sweep.mode must be one of {list(SWEEP_MODES)}, got {mode!r}")

        n_firings = validate_n_firings(data.get("n_firings", 1))
        duration = float(data.get("firing_duration_s", 1.0))
        if duration <= 0.0:
            raise StudyConfigError(
                f"sweep.firing_duration_s must be positive, got {duration}")
        thrusters = tuple(int(t) for t in data.get("thrusters", [1]))
        if not thrusters:
            raise StudyConfigError("sweep.thrusters must not be empty")

        raw_firings = data.get("firings")
        firings: tuple[PrescribedFiringSpec, ...] = ()
        if raw_firings:
            firings = tuple(
                PrescribedFiringSpec.from_mapping(entry, i, thrusters, duration)
                for i, entry in enumerate(raw_firings))
            n_poses = (len(angles) * len(distances)
                       * len(offsets_u) * len(offsets_v))
            expected = n_firings if mode == "per_case" else n_poses * n_firings
            if len(firings) != expected:
                where = ("per case" if mode == "per_case"
                         else f"across {n_poses} poses")
                raise StudyConfigError(
                    f"sweep.n_firings is {n_firings} ({expected} entries "
                    f"{where} in {mode!r} mode) but {len(firings)} explicit "
                    "firings were supplied; the counts must agree exactly")

        return cls(plate_angles_deg=angles, source_distances=distances,
                   n_firings=n_firings, firing_duration_s=duration,
                   thrusters=thrusters, firings=firings, mode=mode,
                   source_offsets_u=offsets_u, source_offsets_v=offsets_v,
                   source_axis_mode=axis_mode)


def validate_n_firings(value: Any) -> int:
    """Validate a requested firing count: a positive integer, nothing else.

    ``n_firings`` is the exact number of Jet Firing History entries a case
    writes, so zero, negative, fractional and non-numeric values are rejected
    rather than coerced. The rule itself lives with the JFH-generation code
    (:func:`pyrpod.rpod.approach_maneuvers.validate_n_firings`) so the
    prescribed and dynamics-driven paths cannot drift apart; only the raised
    error type is specialized here.
    """
    try:
        return _validate_n_firings(value)
    except ValueError as exc:
        raise StudyConfigError(str(exc)) from exc


@dataclass(frozen=True)
class Normalization:
    """Coefficient normalization inputs.

    Coefficients are computed only when the configuration supplies every
    value a given coefficient needs; nothing is invented or defaulted:

    * force coefficients need ``reference_area`` and ``dynamic_pressure``;
    * moment coefficients additionally need ``reference_length``;
    * surface-load coefficients (pressure/shear/heat flux) need
      ``dynamic_pressure`` and, for heat flux, ``reference_heat_flux``.
    """

    reference_area: float | None = None
    reference_length: float | None = None
    dynamic_pressure: float | None = None
    reference_heat_flux: float | None = None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None) -> "Normalization":
        if not data:
            return cls()

        def optional(key: str) -> float | None:
            value = data.get(key)
            if value is None:
                return None
            number = float(value)
            if not np.isfinite(number) or number <= 0.0:
                raise StudyConfigError(
                    f"loads.normalization.{key} must be a positive finite "
                    f"number, got {value!r}")
            return number

        return cls(reference_area=optional("reference_area"),
                   reference_length=optional("reference_length"),
                   dynamic_pressure=optional("dynamic_pressure"),
                   reference_heat_flux=optional("reference_heat_flux"))

    @property
    def has_force_inputs(self) -> bool:
        return (self.reference_area is not None
                and self.dynamic_pressure is not None)

    @property
    def has_moment_inputs(self) -> bool:
        return self.has_force_inputs and self.reference_length is not None

    def to_dict(self) -> dict[str, float | None]:
        return {"reference_area": self.reference_area,
                "reference_length": self.reference_length,
                "dynamic_pressure": self.dynamic_pressure,
                "reference_heat_flux": self.reference_heat_flux}


#: Reference-length mode names for :class:`KnudsenSpec`.
KNUDSEN_REFERENCE_MODES = ("source_distance", "explicit")

#: Definition label implied by each reference-length mode when the
#: configuration does not supply one of its own.
_KNUDSEN_DEFAULT_DEFINITIONS = {
    "source_distance": "lambda_over_source_distance",
    "explicit": "lambda_over_reference_length",
}


@dataclass(frozen=True)
class KnudsenSpec:
    """Derived Knudsen-number METADATA. Never an input to the physics.

    PyRPOD's plume models are collisionless, and this block does not change
    that: no solution, field value or surface load anywhere in the pipeline
    depends on Kn. It exists so an analytical case can be LABELLED with the
    rarefaction regime it is meant to represent, which is what a later,
    entirely separate workflow needs in order to line PyRPOD cases up with
    externally generated DSMC runs.

    ``Kn = mean_free_path_m / reference_length``, with the reference length
    chosen by exactly one of two mutually exclusive modes:

    * ``reference_length: source_distance`` -- the case's own swept source
      distance, so Kn varies across a distance sweep;
    * ``reference_length_m: <number>`` -- a fixed length (a nozzle diameter,
      a panel chord, ...), so Kn is the same for every case.

    The mean free path is always supplied by the configuration. It is never
    inferred from the gas properties in the thruster definition file, because
    a free-molecular model carries no collision rate to infer it from.

    Attributes
    ----------
    mean_free_path_m : float
        Ambient/reference molecular mean free path (m); positive and finite.
    reference_mode : str
        ``'source_distance'`` or ``'explicit'``.
    reference_length_m : float or None
        The fixed reference length, in ``'explicit'`` mode only.
    definition : str
        Free-text label recorded with every result, e.g.
        ``lambda_over_nozzle_diameter``.
    """

    mean_free_path_m: float
    reference_mode: str
    reference_length_m: float | None = None
    definition: str = "lambda_over_source_distance"

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None
                     ) -> "KnudsenSpec | None":
        """Parse the optional ``knudsen`` block; None when it is absent."""
        if not data:
            return None
        if not isinstance(data, Mapping):
            raise StudyConfigError("'knudsen' section must be a mapping")

        if "mean_free_path_m" not in data or data["mean_free_path_m"] is None:
            raise StudyConfigError(
                "knudsen.mean_free_path_m is required; PyRPOD never infers a "
                "mean free path from gas properties (the plume models are "
                "collisionless)")
        try:
            mean_free_path = float(data["mean_free_path_m"])
        except (TypeError, ValueError) as exc:
            raise StudyConfigError(
                "knudsen.mean_free_path_m must be a positive finite number, "
                f"got {data['mean_free_path_m']!r}") from exc
        if not np.isfinite(mean_free_path) or mean_free_path <= 0.0:
            raise StudyConfigError(
                "knudsen.mean_free_path_m must be a positive finite number, "
                f"got {data['mean_free_path_m']!r}")

        symbolic = data.get("reference_length")
        explicit = data.get("reference_length_m")
        if (symbolic is None) == (explicit is None):
            raise StudyConfigError(
                "knudsen requires EXACTLY ONE reference-length mode: either "
                "reference_length: source_distance (the swept distance) or "
                "reference_length_m: <number> (a fixed length); got "
                f"reference_length={symbolic!r}, "
                f"reference_length_m={explicit!r}")

        if symbolic is not None:
            if str(symbolic) != "source_distance":
                raise StudyConfigError(
                    "knudsen.reference_length must be 'source_distance'; for "
                    "any other reference length use reference_length_m: "
                    f"<number>, got {symbolic!r}")
            mode, reference_length = "source_distance", None
        else:
            # The XOR check above guarantees `explicit` is not None here.
            try:
                reference_length = float(explicit)  # type: ignore[arg-type]
            except (TypeError, ValueError) as exc:
                raise StudyConfigError(
                    "knudsen.reference_length_m must be a positive finite "
                    f"number, got {explicit!r}") from exc
            if not np.isfinite(reference_length) or reference_length <= 0.0:
                raise StudyConfigError(
                    "knudsen.reference_length_m must be a positive finite "
                    f"number, got {explicit!r}")
            mode = "explicit"

        definition = data.get("definition")
        return cls(mean_free_path_m=mean_free_path, reference_mode=mode,
                   reference_length_m=reference_length,
                   definition=(str(definition) if definition
                               else _KNUDSEN_DEFAULT_DEFINITIONS[mode]))

    # ------------------------------------------------------------ evaluation
    def reference_length_for(self, source_distance: float) -> float:
        """Reference length used for one case, in this spec's mode."""
        if self.reference_mode == "source_distance":
            distance = float(source_distance)
            if not np.isfinite(distance) or distance <= 0.0:
                raise StudyConfigError(
                    "knudsen.reference_length: source_distance needs a "
                    f"positive finite source distance, got {source_distance!r}")
            return distance
        # 'explicit' mode validates reference_length_m at parse time.
        return float(self.reference_length_m)  # type: ignore[arg-type]

    def knudsen_number(self, source_distance: float) -> float:
        """Derived ``Kn = lambda / L_ref`` for one case. Metadata only."""
        return self.mean_free_path_m / self.reference_length_for(
            source_distance)

    def to_dict(self) -> dict[str, Any]:
        """Plain-data form recorded in the study metadata."""
        return {
            "mean_free_path_m": self.mean_free_path_m,
            "reference_mode": self.reference_mode,
            "reference_length_m": self.reference_length_m,
            "definition": self.definition,
            "role": "derived metadata only; the plume models are "
                    "collisionless and no solution depends on Kn",
        }


@dataclass(frozen=True)
class LoadsSpec:
    """Surface-load integration settings."""

    moment_reference_point: NDArray[np.float64]
    normalization: Normalization

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None) -> "LoadsSpec":
        data = data or {}
        return cls(
            moment_reference_point=_as_vector(
                data.get("moment_reference_point", [0.0, 0.0, 0.0]),
                "loads.moment_reference_point"),
            normalization=Normalization.from_mapping(
                data.get("normalization")),
        )


@dataclass(frozen=True)
class OutputSpec:
    """Output artifact settings (VTK, summary, panel distributions, plots).

    The VTK export stays the PRIMARY full-resolution visualization output and
    is enabled by default. The panel-local surface-distribution CSVs are an
    ADDITIONAL, opt-in export of the same native per-face values in the
    target's own (u, v) coordinates -- convenient for plotting and for a
    later comparison workflow, never a replacement for the VTK files.
    """

    write_vtk: bool = True
    vtk_subdir: str = "vtk"
    summary_csv: str = "case_results.csv"
    summary_metadata: str = "study_metadata.json"
    write_plots: bool = False
    plots_subdir: str = "plots"
    write_surface_distribution: bool = False
    surface_distribution_subdir: str = "distributions"
    write_distribution_plots: bool = False

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None) -> "OutputSpec":
        data = data or {}
        vtk = data.get("vtk") or {}
        summary = data.get("summary") or {}
        plots = data.get("plots") or {}
        distribution = data.get("surface_distribution") or {}
        return cls(
            write_vtk=bool(vtk.get("enabled", True)),
            vtk_subdir=str(vtk.get("subdir", "vtk")),
            summary_csv=str(summary.get("csv", "case_results.csv")),
            summary_metadata=str(summary.get("metadata",
                                             "study_metadata.json")),
            write_plots=bool(plots.get("enabled", False)),
            plots_subdir=str(plots.get("subdir", "plots")),
            write_surface_distribution=bool(distribution.get("enabled",
                                                             False)),
            surface_distribution_subdir=str(distribution.get(
                "subdir", "distributions")),
            # Per-case panel-local pressure maps: only meaningful when the
            # distributions they are drawn from are exported.
            write_distribution_plots=bool(
                plots.get("per_case_distribution", False)),
        )


@dataclass(frozen=True)
class ReferenceSpec:
    """Optional external reference-data location (see mdao.reference_data)."""

    path: str | None = None
    label: str | None = None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None) -> "ReferenceSpec":
        data = data or {}
        path = data.get("path")
        return cls(path=str(path) if path else None,
                   label=str(data["label"]) if data.get("label") else None)


@dataclass(frozen=True)
class StudyConfig:
    """A complete, validated study configuration.

    Attributes
    ----------
    study_name : str
        Identifier recorded with every result row.
    case_dir : str
        Existing PyRPOD case directory (owns ``config.ini`` and the STL / TCD
        / plume assets). Always ends with a path separator, as the rest of
        PyRPOD expects.
    output_dir : str
        Directory the study writes its artifacts to.
    plume_model : str
        The collisionless model that COMPUTES the plume field, one of
        :data:`SUPPORTED_PLUME_MODELS`. Defaults to
        ``'SimplifiedGasKinetics'``.
    knudsen : KnudsenSpec or None
        Optional derived-Knudsen metadata (see :class:`KnudsenSpec`). None
        when the configuration has no ``knudsen`` block, in which case every
        Knudsen field is simply omitted from the results.
    source_path : str
        Path the configuration was read from (configuration provenance).
    """

    study_name: str
    description: str
    case_dir: str
    output_dir: str
    target: TargetSpec
    sweep: SweepSpec
    loads: LoadsSpec
    output: OutputSpec
    reference: ReferenceSpec
    thruster_id: str | None = None
    plume_model: str = SUPPORTED_PLUME_MODEL
    plume_model_parameters: dict[str, Any] = field(default_factory=dict)
    knudsen: KnudsenSpec | None = None
    coordinate_system: str = "case global frame"
    units: dict[str, str] = field(default_factory=lambda: dict(DEFAULT_UNITS))
    source_path: str = ""

    # ------------------------------------------------------------------ load
    @classmethod
    def from_yaml(cls, path: str | os.PathLike[str]) -> "StudyConfig":
        """Read and validate a YAML study configuration."""
        path = os.fspath(path)
        if not os.path.isfile(path):
            raise StudyConfigError(f"study configuration not found: {path!r}")
        with open(path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        if not isinstance(data, Mapping):
            raise StudyConfigError(
                f"study configuration {path!r} must be a YAML mapping")
        return cls.from_mapping(data, source_path=path)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any],
                     source_path: str = "") -> "StudyConfig":
        """Validate an already-parsed configuration mapping."""
        study = data.get("study") or {}
        if not isinstance(study, Mapping):
            raise StudyConfigError("'study' section must be a mapping")

        name = str(_require(study, "name", "study"))
        raw_case_dir = str(_require(study, "case_dir", "study"))
        base_dir = os.path.dirname(os.path.abspath(source_path)) if \
            source_path else os.getcwd()
        case_dir = _resolve_dir(raw_case_dir, base_dir)
        if not os.path.isdir(case_dir):
            raise StudyConfigError(
                f"study.case_dir does not exist: {case_dir!r}")
        if not os.path.isfile(os.path.join(case_dir, "config.ini")):
            raise StudyConfigError(
                f"study.case_dir {case_dir!r} has no config.ini; a study "
                "always runs on top of an existing PyRPOD case")
        case_dir = case_dir.rstrip("\\/") + os.sep

        raw_output_dir = str(study.get(
            "output_dir", os.path.join(case_dir, "results", "studies", name)))
        output_dir = _resolve_dir(raw_output_dir, base_dir)

        plume = data.get("plume_model") or {}
        try:
            model_name = resolve_model_name(
                plume.get("name", SUPPORTED_PLUME_MODEL))
        except PlumeModelError as exc:
            raise StudyConfigError(f"plume_model.name: {exc}") from exc
        model_parameters = dict(plume.get("parameters") or {})

        thruster = data.get("thruster") or {}
        thruster_id = thruster.get("id")

        target_data = data.get("target") or {}
        default_geometry = _default_geometry_id(case_dir)

        metadata = data.get("metadata") or {}
        units = dict(DEFAULT_UNITS)
        units.update({str(k): str(v)
                      for k, v in (metadata.get("units") or {}).items()})

        return cls(
            study_name=name,
            description=str(study.get("description", "")),
            case_dir=case_dir,
            output_dir=output_dir,
            target=TargetSpec.from_mapping(target_data, default_geometry),
            sweep=SweepSpec.from_mapping(data.get("sweep") or {}),
            loads=LoadsSpec.from_mapping(data.get("loads")),
            output=OutputSpec.from_mapping(data.get("output")),
            reference=ReferenceSpec.from_mapping(data.get("reference")),
            thruster_id=str(thruster_id) if thruster_id else None,
            plume_model=model_name,
            plume_model_parameters=model_parameters,
            knudsen=KnudsenSpec.from_mapping(data.get("knudsen")),
            coordinate_system=str(metadata.get("coordinate_system",
                                               "case global frame")),
            units=units,
            source_path=os.path.abspath(source_path) if source_path else "",
        )

    # ------------------------------------------------------------- accessors
    @property
    def n_cases(self) -> int:
        """Number of swept cases: angles x distances x u offsets x v offsets.

        With the default single-element offset lists this is the historical
        angle x distance count; for an offset sweep (which fixes the angle)
        it is ``n_distances * n_u_offsets * n_v_offsets``.
        """
        return len(self.sweep.sweep_poses)

    def with_output_dir(self, output_dir: str) -> "StudyConfig":
        """Copy of this configuration writing to a different directory."""
        return replace(self, output_dir=os.path.abspath(output_dir))

    def provenance(self) -> dict[str, Any]:
        """Configuration provenance recorded in the study metadata."""
        provenance: dict[str, Any] = {
            "study_name": self.study_name,
            "config_path": self.source_path,
            "case_dir": os.path.abspath(self.case_dir),
            "output_dir": os.path.abspath(self.output_dir),
            "plume_model": self.plume_model,
            "plume_model_parameters": dict(self.plume_model_parameters),
            "source_axis_mode": self.sweep.source_axis_mode,
            "coordinate_system": self.coordinate_system,
            "units": dict(self.units),
        }
        if self.knudsen is not None:
            provenance["knudsen"] = self.knudsen.to_dict()
        return provenance


def _resolve_dir(path: str, base_dir: str) -> str:
    """Resolve a configured directory relative to the config file location."""
    if os.path.isabs(path):
        return os.path.abspath(path)
    return os.path.abspath(os.path.join(base_dir, path))


def _default_geometry_id(case_dir: str) -> str:
    """Target geometry id taken from the case's own ``[tv] stl`` entry."""
    import configparser

    config = configparser.ConfigParser()
    config.read(os.path.join(case_dir, "config.ini"))
    if config.has_option("tv", "stl"):
        return str(config["tv"]["stl"])
    return "unknown"
