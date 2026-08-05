"""
Panel-local figures for offset-sweep (ISS-representative) plume studies.

The sibling :mod:`pyrpod.mdao.study_plots` covers the ANGLE sweep: force,
moment and heat flux against approach angle and stand-off distance. This
module covers the OFFSET sweep, in which the plume source is translated
parallel to a flat panel (``sweep.source_axis_mode: parallel_to_normal``)
and the interesting independent variable is the panel-local source offset.

Two families of figure are produced, both optional and both written only
when the study asks for them:

*per-case pressure distribution*
    ``panel_pressure_<case id>.png`` -- the panel-local pressure field of
    one case, drawn from the exported distribution CSV.

*sweep trends*
    ``normal_force_vs_offset_u.png``, ``moment_v_vs_offset_u.png``,
    ``peak_pressure_vs_offset_u.png``, ``cop_u_vs_offset_u.png`` and
    ``normal_force_vs_distance.png`` -- each grouped by stand-off distance
    and plume model, so several models plotted from merged results stay
    visually distinct.

Plotting conventions follow the existing :mod:`pyrpod.mdao.study_plots`:
matplotlib's non-interactive Agg backend is selected on import so a headless
or CI run never opens a window, the same categorical palette is reused, and
no plotting dependency beyond matplotlib is introduced.

A note on the distribution figure
---------------------------------
A flat-plate target is an UNSTRUCTURED triangle mesh, so no structured grid
is fabricated for it. The pressure field is drawn either as a Delaunay
triangulation of the face centroids (filled contours) or, when the face
count is small enough that contours would be misleading, as a face-coloured
scatter of the centroids themselves. Both show the values the strike
pipeline actually computed.
"""

from __future__ import annotations

import csv
import os
from collections import defaultdict
from typing import Iterable, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (backend must be set first)
import matplotlib.tri as mtri    # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402

from pyrpod.mdao.study_plots import PALETTE  # noqa: E402
from pyrpod.mdao.study_results import CaseResult, StudyResults  # noqa: E402

__all__ = [
    "plot_offset_sweep_trends",
    "plot_panel_pressure",
    "plot_panel_pressure_for_case",
]

#: Below this face count a filled-contour plot interpolates more than it
#: reveals, so the faces are drawn individually instead.
MIN_FACES_FOR_CONTOURS = 64

#: Distinct line styles, so two models at the same distance stay separable
#: for a reader who cannot rely on colour alone.
MODEL_STYLES = ("-o", "--s", ":^", "-.d")


# --------------------------------------------------------------------------
# per-case panel-local pressure distribution
# --------------------------------------------------------------------------
def plot_panel_pressure(local_u: Sequence[float], local_v: Sequence[float],
                        pressure: Sequence[float], path: str, *,
                        title: str = "Panel-local pressure distribution",
                        panel_half_u: float | None = None,
                        panel_half_v: float | None = None,
                        centerline_u: float | None = None,
                        centerline_v: float | None = None) -> str:
    """Draw a panel-local pressure field from per-face values.

    Parameters
    ----------
    local_u, local_v : sequence of float
        Panel-local coordinates of the face centroids (m).
    pressure : sequence of float
        Per-face pressure (Pa), in the same order.
    path : str
        Destination PNG path; parent directories are created.
    title : str, optional
        Figure title.
    panel_half_u, panel_half_v : float, optional
        Panel semi-dimensions; when both are given the panel edges are drawn.
        Defaults to the extent of the supplied coordinates.
    centerline_u, centerline_v : float, optional
        Panel-local point where the plume centerline meets the panel, marked
        when supplied.

    Returns
    -------
    str
        The path written.
    """
    u = np.asarray(local_u, dtype=float)
    v = np.asarray(local_v, dtype=float)
    p = np.asarray(pressure, dtype=float)

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    figure, axes = plt.subplots(figsize=(8, 5))

    # A degenerate field (every face identical, or too few faces) would make
    # contouring meaningless, so fall back to drawing the faces themselves.
    use_contours = (u.size >= MIN_FACES_FOR_CONTOURS
                    and float(np.ptp(u)) > 0.0 and float(np.ptp(v)) > 0.0
                    and float(np.ptp(p)) > 0.0)
    mappable: ScalarMappable
    if use_contours:
        triangulation = mtri.Triangulation(u, v)
        mappable = axes.tricontourf(triangulation, p, levels=24,
                                    cmap="viridis")
    else:
        mappable = axes.scatter(u, v, c=p, s=28, cmap="viridis",
                                edgecolors="none")
    bar = figure.colorbar(mappable, ax=axes)
    bar.set_label("pressure (Pa)")

    half_u = panel_half_u if panel_half_u is not None else float(np.max(np.abs(u)))
    half_v = panel_half_v if panel_half_v is not None else float(np.max(np.abs(v)))
    axes.plot([-half_u, half_u, half_u, -half_u, -half_u],
              [-half_v, -half_v, half_v, half_v, -half_v],
              color="0.25", lw=1.2, label="panel edge")

    if centerline_u is not None and centerline_v is not None:
        axes.plot([centerline_u], [centerline_v], marker="x", ms=11, mew=2.2,
                  color="#e8382a", linestyle="none",
                  label="plume centerline")

    axes.set_xlabel("panel-local u (m), longitudinal")
    axes.set_ylabel("panel-local v (m), transverse")
    # On the figure, not the axes: a wide panel is drawn to scale, so an
    # axes title would run under the colorbar.
    figure.suptitle(title, fontsize=9, y=0.98)
    axes.set_aspect("equal", adjustable="box")
    axes.legend(fontsize=8, loc="upper right")
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return path


def plot_panel_pressure_for_case(case: CaseResult, out_dir: str,
                                 panel_half_u: float | None = None,
                                 panel_half_v: float | None = None,
                                 ) -> str | None:
    """Draw one case's pressure distribution from its exported CSV.

    Returns None (rather than raising) when the case carries no distribution
    export, so a study with distributions disabled simply produces no
    per-case figures.
    """
    path = case.surface_distribution_path
    if not path or not os.path.isfile(path):
        return None

    local_u: list[float] = []
    local_v: list[float] = []
    pressure: list[float] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            local_u.append(float(row["local_u"]))
            local_v.append(float(row["local_v"]))
            pressure.append(float(row["pressure"]))
    if not pressure:
        return None

    # In parallel_to_normal mode the plume axis is anti-parallel to the panel
    # normal, so the centerline meets the panel exactly at the source offset.
    centerline = (case.source_offset_u, case.source_offset_v) \
        if case.source_axis_mode == "parallel_to_normal" else (None, None)

    title = (f"Panel-local pressure -- {case.case_id}\n"
             f"{case.model_variant}, L = {case.source_distance:g} m, "
             f"source u = {case.source_offset_u:g} m, "
             f"v = {case.source_offset_v:g} m")
    return plot_panel_pressure(
        local_u, local_v, pressure,
        os.path.join(out_dir, f"panel_pressure_{case.case_id}.png"),
        title=title, panel_half_u=panel_half_u, panel_half_v=panel_half_v,
        centerline_u=centerline[0], centerline_v=centerline[1])


# --------------------------------------------------------------------------
# offset-sweep trends
# --------------------------------------------------------------------------
def _finite(cases: Iterable[CaseResult], quantity: str
            ) -> list[tuple[CaseResult, float]]:
    """Cases whose quantity is present and finite, with that value."""
    selected = []
    for case in cases:
        value = getattr(case, quantity, None)
        if value is None:
            continue
        number = float(value)
        if np.isfinite(number):
            selected.append((case, number))
    return selected


def _grouped_series(cases: Sequence[CaseResult], quantity: str,
                    variable: str, group_by: Sequence[str],
                    ) -> list[tuple[str, list[float], list[float]]]:
    """Series of (label, xs, ys), one per distinct combination of group_by."""
    groups: dict[tuple, list[tuple[float, float]]] = defaultdict(list)
    for case, value in _finite(cases, quantity):
        key = tuple(getattr(case, name) for name in group_by)
        groups[key].append((float(getattr(case, variable)), value))

    series: list[tuple[str, list[float], list[float]]] = []
    for key in sorted(groups, key=lambda k: tuple(str(part) for part in k)):
        points = sorted(groups[key])
        label = ", ".join(
            f"{_label_for(name)}{part:g}" if isinstance(part, (int, float))
            else str(part)
            for name, part in zip(group_by, key))
        series.append((label, [x for x, _ in points], [y for _, y in points]))
    return series


def _label_for(field: str) -> str:
    return {"source_distance": "L = ", "source_offset_u": "u = ",
            "source_offset_v": "v = "}.get(field, f"{field} = ")


def _trend_figure(path: str, title: str, xlabel: str, ylabel: str,
                  series: Sequence[tuple[str, list[float], list[float]]],
                  ) -> str | None:
    """Write one grouped trend figure; None when there is nothing to draw."""
    if not any(xs for _, xs, _ in series):
        return None
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    figure, axes = plt.subplots(figsize=(7, 4.5))
    for index, (label, xs, ys) in enumerate(series):
        axes.plot(xs, ys, MODEL_STYLES[index % len(MODEL_STYLES)],
                  color=PALETTE[index % len(PALETTE)], lw=2, ms=4, label=label)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title, fontsize=10)
    axes.grid(True, color="0.9", lw=0.8)
    axes.legend(fontsize=8)
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return path


def plot_offset_sweep_trends(results: StudyResults, out_dir: str,
                             component: str | None = None) -> list[str]:
    """Write the offset-sweep trend figures; returns the paths written.

    Each figure groups by stand-off distance AND plume model, so results
    merged from a Simplified run and a Collisionless run of the same
    geometry plot as separate, labelled series rather than being averaged
    together.

    Parameters
    ----------
    results : StudyResults
        Executed study results.
    out_dir : str
        Directory to write the figures to (created when missing).
    component : str, optional
        Restrict to one target component; defaults to the first present.
    """
    cases = list(results.cases)
    if not cases:
        return []
    component = component or cases[0].component
    cases = [case for case in cases if case.component == component]
    # One point per pose: keep the first firing of each.
    seen: set[tuple[float, float, float, float, str]] = set()
    unique: list[CaseResult] = []
    for case in cases:
        key = (case.plate_angle_deg, case.source_distance,
               case.source_offset_u, case.source_offset_v, case.plume_model)
        if key in seen:
            continue
        seen.add(key)
        unique.append(case)
    cases = unique
    if not cases:
        return []

    os.makedirs(out_dir, exist_ok=True)
    written: list[str] = []
    group = ("source_distance", "model_variant")

    for quantity, filename, ylabel, title in (
            ("normal_force", "normal_force_vs_offset_u.png",
             "normal force (N), + into panel",
             "Normal force vs source u offset"),
            ("local_moment_v", "moment_v_vs_offset_u.png",
             "local moment about v (N*m)",
             "Panel moment about v vs source u offset"),
            ("max_pressure", "peak_pressure_vs_offset_u.png",
             "peak pressure (Pa)", "Peak pressure vs source u offset"),
            ("center_of_pressure_u", "cop_u_vs_offset_u.png",
             "center of pressure u (m)",
             "Center-of-pressure u vs source u offset")):
        path = _trend_figure(
            os.path.join(out_dir, filename), f"{title} ({component})",
            "source u offset (m)", ylabel,
            _grouped_series(cases, quantity, "source_offset_u", group))
        if path:
            written.append(path)

    # Normal force against stand-off, one series per (u offset, model).
    path = _trend_figure(
        os.path.join(out_dir, "normal_force_vs_distance.png"),
        f"Normal force vs source distance ({component})",
        "source distance L (m)", "normal force (N), + into panel",
        _grouped_series(cases, "normal_force", "source_distance",
                        ("source_offset_u", "model_variant")))
    if path:
        written.append(path)

    # The transverse offset gets its own figures only when it is swept.
    if len({case.source_offset_v for case in cases}) > 1:
        for quantity, filename, ylabel, title in (
                ("normal_force", "normal_force_vs_offset_v.png",
                 "normal force (N), + into panel",
                 "Normal force vs source v offset"),
                ("local_moment_u", "moment_u_vs_offset_v.png",
                 "local moment about u (N*m)",
                 "Panel moment about u vs source v offset"),
                ("center_of_pressure_v", "cop_v_vs_offset_v.png",
                 "center of pressure v (m)",
                 "Center-of-pressure v vs source v offset")):
            path = _trend_figure(
                os.path.join(out_dir, filename), f"{title} ({component})",
                "source v offset (m)", ylabel,
                _grouped_series(cases, quantity, "source_offset_v", group))
            if path:
                written.append(path)

    return written
