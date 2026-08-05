"""INDEPENDENT Cai 2016 analytical reference for the ISS-panel case.

This is a **reference generator, not a plume-model backend.** It is not
imported by PyRPOD, is not reachable from `PlumeStrikeCalculator`, and takes
no part in any study run. Its only job is to evaluate the *exact* Cai 2016
surface solution (`pyrpod/plume/CaiImpingement2016.py`, Eqs. 9-14 by
Gauss-Legendre quadrature over the nozzle exit disk) at the **same face
centroids** a study used, and export it in the **same panel-local
distribution schema**, so the two can be diffed column for column.

Why it is worth having
----------------------
The production pipeline reaches the wall through a chain -- plume field ->
`LocalFieldState` -> Maxwellian (Shen) wall formulas -- evaluated per face.
The Cai 2016 surface solution integrates the incident *and* re-emitted
molecular fluxes at the wall directly. Agreement between them is an
independent check on the chain; disagreement localizes to it. That check
already exists for integrated loads in `tests/mdao/mdao_integration_test_02.py`;
this script provides it per face.

Scope and limits
----------------
* **Normal incidence only.** The exported geometry must be a flat panel with
  the plume axis anti-parallel to its normal, i.e. a study running
  `source_axis_mode: parallel_to_normal`. That is Cai's `alpha_0 = 90 deg`.
  Anything else is refused rather than silently mapped.
* **Diffuse coefficients** (`Cp_d`, `Cf1_d`, `Cf2_d`, `Cq_d`) are exported;
  the case is fully diffuse (`sigma = 1`).
* **No DSMC.** This is an analytical reference, nothing more.

Usage
-----
    # Evaluate at the committed panel mesh, for the baseline pose
    python cai2016_reference.py --distance 4.0

    # Match a swept case: source offset 5.5 m along u
    python cai2016_reference.py --distance 4.0 --offset-u 5.5

    # Compare against a distribution a study already exported
    python cai2016_reference.py --distance 4.0 \
        --compare results/studies/baseline_simplified/cases/<case>/distributions/<case>_panel_firing001.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parents[2]))

from pyrpod.mdao.surface_distribution import (  # noqa: E402
    DISTRIBUTION_COLUMNS,
    write_surface_distribution,
)
from pyrpod.plume import CaiImpingement2016 as cai  # noqa: E402

# Case conditions, matching config.ini and tcd/tdf.csv (Cai 2016 Section 4).
S_0 = 2.0                       # exit speed ratio
EPS = 1.5                       # Tw / T0 = 300 / 200
R_0 = 0.5                       # nozzle exit radius (m), D = 1 m
Q_DYN = 1.1044652197738332      # n0*m*U0^2/2 (Pa)
Q_HEAT = 637.3520362956127      # n0*m*U0^3/2 (W/m^2)

DEFAULT_STL = CASE_DIR / "stl" / "iss_panel.stl"


def reference_distribution(centroids: np.ndarray, areas: np.ndarray,
                           distance: float, offset_u: float = 0.0,
                           offset_v: float = 0.0) -> list[dict[str, float]]:
    """Exact Cai 2016 surface loads at panel face centroids.

    Frame mapping. The panel lies in the global X-Y plane with its normal
    along +Z; the source stands off at ``+distance * n`` and is translated to
    ``(offset_u, offset_v)``, with its axis along ``-n``. Cai's nozzle frame
    puts the exit at the origin with the jet along ``+X`` and, for
    ``alpha_0 = 90 deg``, the plate in the ``Y-Z`` plane at ``X = L``:

        X_paper = distance                     (along the jet axis)
        Y_paper = v - offset_v                 (transverse, Cai's s)
        Z_paper = u - offset_u                 (longitudinal, Cai's tau)

    The longitudinal panel axis is mapped to Cai's ``tau`` so that ``Cf1_d``,
    the shear along ``tau``, is the longitudinal shear. At normal incidence
    the solution is axisymmetric, so the choice affects only which shear
    component is which, never the exported magnitudes.

    Returns
    -------
    list of dict
        Rows in the :data:`DISTRIBUTION_COLUMNS` schema. ``strike_count`` is
        1.0 for every face carrying a non-zero load: this analytical solution
        has no wedge/radius gating, so it is a coverage flag, not a tally.
    """
    centroids = np.asarray(centroids, dtype=float)
    local_u = centroids[:, 0]
    local_v = centroids[:, 1]

    X = np.full(local_u.shape, float(distance))
    Y = local_v - float(offset_v)
    Z = local_u - float(offset_u)

    field = cai.surface_coefficients(X, Y, Z, S_0, np.pi / 2.0, EPS, R_0)

    pressure = np.asarray(field["Cp_d"], dtype=float) * Q_DYN
    shear = np.hypot(np.asarray(field["Cf1_d"], dtype=float),
                     np.asarray(field["Cf2_d"], dtype=float)) * Q_DYN
    heat_flux = np.asarray(field["Cq_d"], dtype=float) * Q_HEAT
    loaded = (pressure != 0.0) | (shear != 0.0) | (heat_flux != 0.0)

    return [
        {
            "face_index": int(i),
            "centroid_x": float(centroids[i, 0]),
            "centroid_y": float(centroids[i, 1]),
            "centroid_z": float(centroids[i, 2]),
            "local_u": float(local_u[i]),
            "local_v": float(local_v[i]),
            "area": float(areas[i]),
            "pressure": float(pressure[i]),
            "shear_stress": float(shear[i]),
            "heat_flux": float(heat_flux[i]),
            "strike_count": 1.0 if loaded[i] else 0.0,
        }
        for i in range(len(centroids))
    ]


def compare_with(rows: list[dict[str, float]], path: str) -> None:
    """Report per-face agreement against a study-exported distribution CSV."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        study = list(csv.DictReader(handle))
    if len(study) != len(rows):
        raise SystemExit(
            f"face-count mismatch: reference has {len(rows)}, "
            f"{path} has {len(study)}; both must come from the same mesh")

    print(f"\ncompared against {path}")
    print(f"{'quantity':14s} {'ref max':>12s} {'study max':>12s} "
          f"{'max |diff|':>12s} {'rel. of peak':>13s}")
    for column in ("pressure", "shear_stress", "heat_flux"):
        reference = np.array([row[column] for row in rows], dtype=float)
        candidate = np.array([float(row[column]) for row in study],
                             dtype=float)
        difference = np.abs(reference - candidate)
        peak = float(np.max(np.abs(reference)))
        relative = float(np.max(difference) / peak) if peak > 0.0 else float("nan")
        print(f"{column:14s} {peak:12.6g} {float(np.max(candidate)):12.6g} "
              f"{float(np.max(difference)):12.6g} {relative:13.4%}")
    print("\nA gap here is the Maxwellian wall chain, not a bug: PyRPOD "
          "reaches the wall through plume field -> LocalFieldState -> Shen "
          "formulas, while this reference integrates the wall fluxes "
          "directly.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stl", type=Path, default=DEFAULT_STL,
                        help="panel STL to evaluate at")
    parser.add_argument("--distance", type=float, required=True,
                        help="source stand-off L along the panel normal (m)")
    parser.add_argument("--offset-u", type=float, default=0.0,
                        help="source offset along the longitudinal axis (m)")
    parser.add_argument("--offset-v", type=float, default=0.0,
                        help="source offset along the transverse axis (m)")
    parser.add_argument("--out", type=Path, default=None,
                        help="output CSV (default: cai2016_reference_*.csv "
                             "beside this script)")
    parser.add_argument("--compare", type=Path, default=None,
                        help="a study-exported distribution CSV to diff against")
    args = parser.parse_args(argv)

    if args.distance <= 0.0:
        raise SystemExit(f"--distance must be positive, got {args.distance}")

    from stl import mesh

    panel = mesh.Mesh.from_file(str(args.stl))
    normals = panel.get_unit_normals()
    if not np.allclose(normals, [0.0, 0.0, 1.0], atol=1e-6):
        raise SystemExit(
            f"{args.stl} is not a flat panel with a +Z normal; this "
            "reference generator implements normal incidence "
            "(Cai alpha_0 = 90 deg) only")

    centroids = panel.vectors.mean(axis=1)
    v0, v1, v2 = panel.vectors[:, 0], panel.vectors[:, 1], panel.vectors[:, 2]
    areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)

    rows = reference_distribution(centroids, areas, args.distance,
                                  args.offset_u, args.offset_v)

    out = args.out or (CASE_DIR / f"cai2016_reference_L{args.distance:g}"
                       f"_u{args.offset_u:g}_v{args.offset_v:g}.csv")
    write_surface_distribution(out, rows, {
        "source": "INDEPENDENT analytical reference; NOT a PyRPOD plume "
                  "model and not used by any study run",
        "solution": "Cai 2016 (Aerospace 3(4):43) Eqs. 9-14, exact "
                    "quadrature over the nozzle exit disk",
        "generator": "case/plume/iss_panel_thesis/cai2016_reference.py",
        "geometry": str(args.stl),
        "alpha_0_deg": 90.0,
        "source_distance": float(args.distance),
        "source_offset_u": float(args.offset_u),
        "source_offset_v": float(args.offset_v),
        "speed_ratio_S0": S_0,
        "temperature_ratio_eps": EPS,
        "nozzle_radius_m": R_0,
        "dynamic_pressure_Pa": Q_DYN,
        "reference_heat_flux_W_m2": Q_HEAT,
        "accommodation": "fully diffuse (sigma = 1); diffuse coefficients "
                         "Cp_d, Cf1_d, Cf2_d, Cq_d",
        "columns_note": f"schema matches {list(DISTRIBUTION_COLUMNS)}",
    })
    print(f"wrote {out} ({len(rows)} faces)")

    if args.compare:
        compare_with(rows, str(args.compare))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
