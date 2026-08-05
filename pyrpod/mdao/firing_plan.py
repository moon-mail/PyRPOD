"""
Prescribed Jet Firing History generation for plume/target validation studies.

Validation sweeps do not fly a trajectory: the firing poses are PRESCRIBED
(placed by the engineer, or generated from a swept approach angle and source
distance about a stationary target) and the Jet Firing History is written
from them directly. This module owns that generation, so the dynamics-based
approach-maneuver workflow in ``pyrpod.rpod.approach_maneuvers`` stays
untouched.

``n_firings`` semantics
-----------------------
``n_firings`` is the EXACT number of entries written to the JFH. The count is
validated before anything is generated (:func:`validate_n_firings`), the
generated sequence is asserted to match it (:func:`build_case_firings`), and
an explicitly prescribed firing list whose length disagrees is an error --
never a silent truncation or extension.

Pose conventions
----------------
Two axis modes are available; a study picks one with
``sweep.source_axis_mode`` (see
:data:`pyrpod.mdao.study_config.SOURCE_AXIS_MODES`). Both use the target
reference point ``C``, its outward normal ``n_hat`` (pointing toward the
plume-source side) and its in-plane tangent ``u_hat``, together with the
transverse axis ``v_hat = n_hat x u_hat``.

``aim_at_reference`` (default, unchanged)
    The source sits on an arc about ``C`` and always aims back at it:

        d_hat(alpha) = cos(alpha) * n_hat + sin(alpha) * u_hat
        source position = C + L * d_hat(alpha)
        thruster axis   = -d_hat(alpha)            (aimed at C)

    ``alpha = 0`` is head-on. This reproduces the pose convention of the
    committed sweep-JFH generators
    (``case/plume/plume_flat_plate_sweep/jfh/generate_sweep_jfh.py``)
    exactly, including the binormal choice ``cross(n_hat, u_hat)`` for the
    DCM's second column.

``parallel_to_normal`` (ISS-panel studies)
    The source is TRANSLATED parallel to the surface and its axis stays
    fixed, so the plume centerline strikes the panel at the requested
    panel-local offset instead of always at ``C``:

        source position = C + L * n_hat + u_off * u_hat + v_off * v_hat
        thruster axis   = -n_hat

    This is deliberately NOT the same experiment as moving the source while
    continuously re-aiming it at the panel center. At zero offsets the two
    modes coincide (``aim_at_reference`` at ``alpha = 0``), so the new mode
    is a strict extension of the old convention rather than a redefinition.

In both modes the JFH DCM is built with the thruster axis as its first
COLUMN, which is what the strike pipeline reads as the plume normal
(``dcm.T`` rows, with an identity thruster DCM).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import NDArray

from pyrpod.mdao.study_config import (
    PrescribedFiringSpec,
    StudyConfigError,
    SweepPose,
    SweepSpec,
    TargetSpec,
    validate_n_firings,
)

__all__ = [
    "Firing",
    "build_case_firings",
    "build_sweep_firings",
    "pose_for",
    "pose_for_sweep_pose",
    "translated_pose_for",
    "validate_n_firings",
    "write_jfh_file",
]


@dataclass(frozen=True)
class Firing:
    """One JFH entry: pose, active thrusters, firing duration and start time.

    Attributes
    ----------
    position : np.ndarray
        Visiting-vehicle (plume-source) position in the case's global frame.
    dcm : np.ndarray
        3x3 direction cosine matrix; its first column is the thruster axis.
    thrusters : tuple of int
        JFH thruster indices active during this firing (1-based).
    duration_s : float
        Firing time; the strike pipeline multiplies heat flux by this to get
        the per-firing heat-flux load.
    start_time_s : float
        Elapsed time at the start of this firing (informational; the strike
        calculation does not read it).
    plate_angle_deg, source_distance : float or None
        The swept parameters this firing realizes. Set when the firing came
        from a sweep pose, so a history spanning many poses stays keyed to
        the sweep grid; None for firings whose pose was prescribed outright
        without a sweep parameterization.
    pose_index : int or None
        Index of the firing's pose in the sweep's execution order.
    source_offset_u, source_offset_v : float
        Panel-local translation of the plume source along the target's
        longitudinal and transverse axes. Zero unless the pose came from an
        offset sweep.
    source_axis_mode : str
        Which pose convention built this firing (see the module docstring).
    """

    position: NDArray[np.float64]
    dcm: NDArray[np.float64]
    thrusters: tuple[int, ...]
    duration_s: float = 1.0
    start_time_s: float = 0.0
    plate_angle_deg: float | None = None
    source_distance: float | None = None
    pose_index: int | None = None
    source_offset_u: float = 0.0
    source_offset_v: float = 0.0
    source_axis_mode: str = "aim_at_reference"


def pose_for(alpha_deg: float, distance: float,
             reference_point: Sequence[float] | NDArray[np.float64],
             normal: Sequence[float] | NDArray[np.float64],
             tangent: Sequence[float] | NDArray[np.float64],
             ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Plume-source pose for one (approach angle, distance) combination.

    Parameters
    ----------
    alpha_deg : float
        Approach angle in degrees; 0 is head-on along ``normal``, positive
        angles rotate the source toward ``tangent``.
    distance : float
        Distance from ``reference_point`` to the plume source.
    reference_point, normal, tangent : array-like
        Target geometry axes (see :class:`pyrpod.mdao.study_config.TargetSpec`).

    Returns
    -------
    (np.ndarray, np.ndarray)
        The source position and the 3x3 DCM whose first column is the
        thruster axis (aimed back at ``reference_point``).
    """
    center = np.asarray(reference_point, dtype=float)
    n_hat = np.asarray(normal, dtype=float)
    t_hat = np.asarray(tangent, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    t_hat = t_hat / np.linalg.norm(t_hat)

    alpha = np.deg2rad(float(alpha_deg))
    d_hat = np.cos(alpha) * n_hat + np.sin(alpha) * t_hat
    position = center + float(distance) * d_hat
    axis = -d_hat                       # aimed at the target reference point

    # Right-handed triad completed with the sweep plane's binormal, matching
    # the committed sweep-JFH generators (which hard-code the global Y axis
    # for their X-Z sweep plane -- the same vector this expression yields).
    binormal = np.cross(n_hat, t_hat)
    binormal = binormal / np.linalg.norm(binormal)
    dcm = np.column_stack([axis, binormal, np.cross(axis, binormal)])
    return position, dcm


def translated_pose_for(distance: float, offset_u: float, offset_v: float,
                        reference_point: Sequence[float] | NDArray[np.float64],
                        normal: Sequence[float] | NDArray[np.float64],
                        tangent: Sequence[float] | NDArray[np.float64],
                        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Plume-source pose TRANSLATED parallel to the target surface.

    The source stands off by ``distance`` along the target normal and is
    then slid along the surface-local axes; its axis stays anti-parallel to
    the normal, so the plume centerline intersects the surface at the
    panel-local point ``(offset_u, offset_v)`` rather than at the reference
    point:

        position = C + L * n_hat + offset_u * u_hat + offset_v * v_hat
        axis     = -n_hat

    Parameters
    ----------
    distance : float
        Stand-off distance L from the surface along ``normal``.
    offset_u, offset_v : float
        Panel-local translations along the longitudinal (``tangent``) and
        transverse (``normal x tangent``) axes.
    reference_point, normal, tangent : array-like
        Target geometry axes (see
        :meth:`pyrpod.mdao.study_config.TargetSpec.local_basis`).

    Returns
    -------
    (np.ndarray, np.ndarray)
        The source position and the 3x3 DCM whose FIRST COLUMN is the
        thruster axis. The triad is the same one
        :func:`pose_for` builds at ``alpha = 0``, so the two modes agree
        exactly when both offsets are zero.
    """
    center = np.asarray(reference_point, dtype=float)
    n_hat = np.asarray(normal, dtype=float)
    u_hat = np.asarray(tangent, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    u_hat = u_hat / np.linalg.norm(u_hat)
    v_hat = np.cross(n_hat, u_hat)
    v_hat = v_hat / np.linalg.norm(v_hat)

    position = (center + float(distance) * n_hat
                + float(offset_u) * u_hat + float(offset_v) * v_hat)
    axis = -n_hat                       # fixed: parallel to -n at every offset

    # Same column convention as pose_for: [axis, binormal, axis x binormal],
    # with the binormal cross(n_hat, u_hat) = v_hat.
    dcm = np.column_stack([axis, v_hat, np.cross(axis, v_hat)])
    return position, dcm


def pose_for_sweep_pose(pose: SweepPose, target: TargetSpec,
                        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Plume-source pose for one :class:`SweepPose`, in its own axis mode.

    The single place the axis mode is turned into geometry, so the per-case
    and single-history engines cannot drift apart.

    Raises
    ------
    StudyConfigError
        If the pose carries an axis mode this module does not implement.
    """
    if pose.axis_mode == "parallel_to_normal":
        return translated_pose_for(pose.source_distance, pose.source_offset_u,
                                   pose.source_offset_v,
                                   target.reference_point, target.normal,
                                   target.tangent)
    if pose.axis_mode == "aim_at_reference":
        return pose_for(pose.plate_angle_deg, pose.source_distance,
                        target.reference_point, target.normal, target.tangent)
    raise StudyConfigError(
        f"unknown sweep.source_axis_mode {pose.axis_mode!r}")


def build_case_firings(sweep: SweepSpec, target: TargetSpec,
                       alpha_deg: float, distance: float,
                       pose_index: int = 0,
                       start_time_s: float = 0.0,
                       pose: SweepPose | None = None) -> list[Firing]:
    """Build EXACTLY ``sweep.n_firings`` firings for one sweep pose.

    When the configuration prescribes firings explicitly they are used
    verbatim (their count is validated against ``n_firings`` when the
    configuration is parsed). Otherwise the swept pose is generated in the
    sweep's own axis mode and repeated for ``n_firings`` successive firing
    intervals -- one JFH entry per requested firing.

    ``pose_index`` and ``start_time_s`` place this pose inside a longer
    sequence; they matter only when many poses share one history (see
    :func:`build_sweep_firings`).

    Parameters
    ----------
    pose : SweepPose, optional
        The fully parameterized pose to realize, including its panel-local
        offsets. When omitted, one is built from ``alpha_deg`` / ``distance``
        with zero offsets in the sweep's axis mode, which is exactly the
        historical behavior. When supplied it is authoritative and
        ``alpha_deg`` / ``distance`` are ignored.

    Raises
    ------
    StudyConfigError
        If the generated sequence length would differ from ``n_firings``.
    """
    n_firings = validate_n_firings(sweep.n_firings)
    if pose is None:
        pose = SweepPose(plate_angle_deg=float(alpha_deg),
                         source_distance=float(distance),
                         axis_mode=sweep.source_axis_mode)

    if sweep.firings:
        specs = sweep.firings
        if sweep.mode == "single_jfh":
            # One shared history: each pose takes its own contiguous slice of
            # the prescribed sequence, in sweep order.
            offset = pose_index * n_firings
            specs = sweep.firings[offset:offset + n_firings]
        firings = [
            _from_spec(spec, index, pose, pose_index, start_time_s)
            for index, spec in enumerate(specs)
        ]
    else:
        position, dcm = pose_for_sweep_pose(pose, target)
        firings = [
            Firing(position=position, dcm=dcm, thrusters=sweep.thrusters,
                   duration_s=sweep.firing_duration_s,
                   start_time_s=start_time_s
                   + index * sweep.firing_duration_s,
                   plate_angle_deg=pose.plate_angle_deg,
                   source_distance=pose.source_distance,
                   pose_index=pose_index,
                   source_offset_u=pose.source_offset_u,
                   source_offset_v=pose.source_offset_v,
                   source_axis_mode=pose.axis_mode)
            for index in range(n_firings)
        ]

    if len(firings) != n_firings:
        raise StudyConfigError(
            f"requested n_firings={n_firings} but built {len(firings)} "
            "firings; the two must agree exactly")
    return firings


def build_sweep_firings(sweep: SweepSpec,
                        target: TargetSpec) -> list[Firing]:
    """Build the WHOLE sweep as one firing sequence.

    Every pose of ``sweep.sweep_poses`` contributes ``sweep.n_firings``
    entries, in execution order, with firing times running continuously
    across the sequence. The result is the single Jet Firing History of a
    ``single_jfh`` study, and its length is exactly ``sweep.total_firings``.

    Raises
    ------
    StudyConfigError
        If the assembled sequence length would differ from that product.
    """
    firings: list[Firing] = []
    elapsed = 0.0
    for pose_index, pose in enumerate(sweep.sweep_poses):
        pose_firings = build_case_firings(
            sweep, target, pose.plate_angle_deg, pose.source_distance,
            pose_index=pose_index, start_time_s=elapsed, pose=pose)
        firings.extend(pose_firings)
        elapsed += sum(firing.duration_s for firing in pose_firings)

    if len(firings) != sweep.total_firings:
        raise StudyConfigError(
            f"sweep of {len(sweep.sweep_poses)} poses x n_firings="
            f"{sweep.n_firings} must produce {sweep.total_firings} JFH "
            f"entries, built {len(firings)}")
    return firings


def _from_spec(spec: PrescribedFiringSpec, index: int,
               pose: SweepPose | None = None,
               pose_index: int | None = None,
               start_time_s: float = 0.0) -> Firing:
    return Firing(position=np.asarray(spec.position, dtype=float),
                  dcm=np.asarray(spec.dcm, dtype=float),
                  thrusters=spec.thrusters,
                  duration_s=spec.duration_s,
                  start_time_s=start_time_s + index * spec.duration_s,
                  plate_angle_deg=(None if pose is None
                                   else pose.plate_angle_deg),
                  source_distance=(None if pose is None
                                   else pose.source_distance),
                  pose_index=pose_index,
                  source_offset_u=(0.0 if pose is None
                                   else pose.source_offset_u),
                  source_offset_v=(0.0 if pose is None
                                   else pose.source_offset_v),
                  source_axis_mode=("aim_at_reference" if pose is None
                                    else pose.axis_mode))


def write_jfh_file(path: str | os.PathLike[str],
                   firings: Iterable[Firing]) -> int:
    """Write firings to a JFH file readable by ``JetFiringHistory.read_jfh``.

    The emitted format matches the committed case generators byte-for-byte in
    structure (header line with the firing count, an unused second line, then
    one row per firing: index, dt, t, unused column, nine DCM values, three
    position values, uncertainty factor, thruster count, thruster indices).

    Returns
    -------
    int
        The number of entries written, which is the number of firings the
        JFH will report.
    """
    firings = list(firings)
    if not firings:
        raise StudyConfigError(
            "refusing to write a JFH with zero firings; n_firings must be a "
            "positive integer")

    path = os.fspath(path)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    lines = [f"offseted    {len(firings)}       0",
             "       0.000       0.000       0.000"]
    for index, firing in enumerate(firings, start=1):
        dcm = np.asarray(firing.dcm, dtype=float).reshape(3, 3)
        dcm_str = " ".join(f"{value:.6e}" for value in dcm.ravel())
        xyz_str = " ".join(f"{value:.9g}"
                           for value in np.asarray(firing.position,
                                                   dtype=float).ravel())
        thruster_str = " ".join(str(int(t)) for t in firing.thrusters)
        lines.append(
            f"      {index} {firing.start_time_s:g} {firing.duration_s:g} 0 "
            f"{dcm_str} {xyz_str} 1 {len(firing.thrusters)} {thruster_str}")

    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(firings)
