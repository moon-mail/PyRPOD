"""
Panel-local surface-distribution export for plume/target validation studies.

The per-firing VTK files remain the PRIMARY full-resolution visualization
output of a study; nothing here replaces or reformats them. This module adds
a second, flat view of exactly the same numbers: one CSV per case and
component holding every face of that component with its panel-local ``(u, v)``
coordinates alongside the native pressure, shear stress, heat flux and strike
count.

Why a second export
-------------------
A ``.vtu`` file is the right artifact for ParaView and the wrong one for a
plotting script, a spreadsheet, or a later comparison workflow that needs to
line PyRPOD faces up with an externally generated dataset. The CSV carries
the same values in the target's own coordinates, so a panel study can be
plotted and compared without a VTK reader.

What it is NOT
--------------
* No interpolation, smoothing, resampling or structured-grid projection.
  Every row is one native mesh face, with the value the strike pipeline
  computed for it. A flat-plate mesh is unstructured triangles and is
  exported as unstructured triangles.
* No common-grid projection onto any external mesh, and nothing DSMC-aware.
  Producing a shared grid is a separate workflow's job.

Units and provenance
--------------------
The CSV holds numbers only. A sidecar JSON written beside it records the
units of every column, the panel basis the coordinates were taken on, the
case's pose, the plume model and the derived Knudsen metadata, so a
distribution file is self-describing without the study metadata document.
"""

from __future__ import annotations

import csv
import json
import os
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from pyrpod.mdao.surface_loads import panel_local_coordinates

__all__ = [
    "DISTRIBUTION_COLUMNS",
    "DISTRIBUTION_UNITS",
    "distribution_rows",
    "write_surface_distribution",
]

#: Column order of an exported distribution CSV. Flat and stable: a reader
#: may rely on these names existing, in this order.
DISTRIBUTION_COLUMNS: tuple[str, ...] = (
    "face_index",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "local_u",
    "local_v",
    "area",
    "pressure",
    "shear_stress",
    "heat_flux",
    "strike_count",
)

#: SI units of every exported column, recorded in the sidecar JSON.
DISTRIBUTION_UNITS: dict[str, str] = {
    "face_index": "-",
    "centroid_x": "m",
    "centroid_y": "m",
    "centroid_z": "m",
    "local_u": "m",
    "local_v": "m",
    "area": "m^2",
    "pressure": "Pa",
    "shear_stress": "Pa",
    "heat_flux": "W/m^2",
    "strike_count": "-",
}


def distribution_rows(
    face_indices: Sequence[int] | NDArray[np.int64],
    centroids: NDArray[np.float64],
    areas: NDArray[np.float64],
    pressures: NDArray[np.float64],
    shear_stresses: NDArray[np.float64],
    heat_fluxes: NDArray[np.float64],
    strikes: NDArray[np.float64] | None,
    *,
    reference_point: Sequence[float] | NDArray[np.float64],
    u_hat: Sequence[float] | NDArray[np.float64],
    v_hat: Sequence[float] | NDArray[np.float64],
) -> list[dict[str, float]]:
    """Build one distribution row per face of a component.

    The arrays are the FULL target mesh's per-face fields; ``face_indices``
    selects the component's faces and is preserved in the ``face_index``
    column, so a row can always be traced back to the mesh and to the VTK
    file. Values are copied through unchanged.

    Parameters
    ----------
    face_indices : array-like of int
        Component face indices into the full-mesh arrays.
    centroids, areas : np.ndarray
        Full-mesh face centroids (N, 3) and areas (N,).
    pressures, shear_stresses, heat_fluxes : np.ndarray
        Full-mesh per-face fields for one firing (Pa, Pa, W/m^2).
    strikes : np.ndarray or None
        Full-mesh per-face strike counts; zeros are recorded when None.
    reference_point, u_hat, v_hat : array-like
        Panel-local origin and in-surface axes (see
        :meth:`pyrpod.mdao.study_config.TargetSpec.local_basis`).

    Returns
    -------
    list of dict
        One dictionary per face, keyed by :data:`DISTRIBUTION_COLUMNS`.
    """
    indices = np.asarray(face_indices, dtype=np.int64)
    centroid = np.asarray(centroids, dtype=float)[indices]
    area = np.asarray(areas, dtype=float)[indices]
    pressure = np.asarray(pressures, dtype=float)[indices]
    shear = np.asarray(shear_stresses, dtype=float)[indices]
    heat_flux = np.asarray(heat_fluxes, dtype=float)[indices]
    strike = (np.zeros(indices.size) if strikes is None
              else np.asarray(strikes, dtype=float)[indices])

    local_u, local_v = panel_local_coordinates(centroid, reference_point,
                                               u_hat, v_hat)

    return [
        {
            "face_index": int(indices[i]),
            "centroid_x": float(centroid[i, 0]),
            "centroid_y": float(centroid[i, 1]),
            "centroid_z": float(centroid[i, 2]),
            "local_u": float(local_u[i]),
            "local_v": float(local_v[i]),
            "area": float(area[i]),
            "pressure": float(pressure[i]),
            "shear_stress": float(shear[i]),
            "heat_flux": float(heat_flux[i]),
            "strike_count": float(strike[i]),
        }
        for i in range(indices.size)
    ]


def write_surface_distribution(path: str | os.PathLike[str],
                               rows: Sequence[Mapping[str, Any]],
                               metadata: Mapping[str, Any] | None = None,
                               ) -> str:
    """Write one component's distribution CSV and its sidecar JSON.

    Parameters
    ----------
    path : str or path-like
        Destination CSV path; parent directories are created.
    rows : sequence of mapping
        Rows from :func:`distribution_rows`.
    metadata : mapping, optional
        Case metadata (pose, plume model, Knudsen, panel basis, ...) written
        to ``<path stem>.meta.json`` together with the column units. When
        omitted only the units and column order are recorded.

    Returns
    -------
    str
        The CSV path written.

    Raises
    ------
    ValueError
        If ``rows`` is empty -- an empty distribution is a bug in the caller,
        not a valid artifact.
    """
    if not rows:
        raise ValueError(
            "refusing to write an empty surface distribution; the component "
            "selects no faces")

    path = os.fspath(path)
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(DISTRIBUTION_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)

    sidecar = f"{os.path.splitext(path)[0]}.meta.json"
    document: dict[str, Any] = {
        "schema": "pyrpod.surface_distribution/1",
        "columns": list(DISTRIBUTION_COLUMNS),
        "units": dict(DISTRIBUTION_UNITS),
        "n_faces": len(rows),
        "interpolation": "none; every row is one native mesh face",
    }
    if metadata:
        document.update(dict(metadata))
    with open(sidecar, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(document, handle, indent=2, sort_keys=False)
        handle.write("\n")
    return path
