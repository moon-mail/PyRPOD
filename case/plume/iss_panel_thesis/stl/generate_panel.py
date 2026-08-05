"""Generate the ISS-representative solar-panel target mesh.

Builds a single-sided, uniformly triangulated rectangular plate standing in
for one ISS solar-array wing:

    length  22 m   along the panel-local u (longitudinal) axis -> global X
    width   12 m   along the panel-local v (transverse)  axis -> global Y
    normal  +Z,    geometric center at the origin

The plume source stands off on the +Z side and, in the study's
``parallel_to_normal`` axis mode, is translated parallel to the panel, so
every face normal must point back toward the source (+Z). The strike
pipeline's facing test (``surface_dot_plume < 0``) requires exactly that;
the script asserts it before saving.

Meshing
-------
The quad grid comes from `surfmesh <https://github.com/plume-kit/surfmesh>`_
(``quad_faces_from_edges`` + ``convert_2d_face_to_3d``), and each quad is
split into the two triangles an STL needs. surfmesh is a build-time
dependency of THIS SCRIPT only: the generated STL is committed, so neither
PyRPOD nor the automated tests import it.

Resolution is parametrized. The committed default is 44 x 24 quads
(0.5 m elements, 2112 triangles) -- fine enough to resolve the impingement
footprint of a plume standing off a few metres, coarse enough that a full
study runs locally in seconds to a couple of minutes. Research-grade runs
should raise ``--n-u`` / ``--n-v`` and regenerate.

Run from this directory:  python generate_panel.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from stl import mesh

#: ISS-representative solar-array wing dimensions (m).
PANEL_LENGTH_U = 22.0
PANEL_WIDTH_V = 12.0

#: Committed mesh resolution: quads per axis (triangles = 2 * n_u * n_v).
DEFAULT_N_U = 44
DEFAULT_N_V = 24


def build_panel_mesh(length_u: float = PANEL_LENGTH_U,
                     width_v: float = PANEL_WIDTH_V,
                     n_u: int = DEFAULT_N_U, n_v: int = DEFAULT_N_V,
                     center: tuple[float, float, float] = (0.0, 0.0, 0.0),
                     ) -> mesh.Mesh:
    """Return a numpy-stl Mesh of the flat rectangular panel.

    Parameters
    ----------
    length_u : float
        Panel length along the panel-local u axis (global X), in m.
    width_v : float
        Panel width along the panel-local v axis (global Y), in m.
    n_u, n_v : int
        Quad divisions per axis; the mesh holds ``2 * n_u * n_v`` triangles.
    center : tuple of float
        Panel geometric center in the global frame.

    Returns
    -------
    stl.mesh.Mesh
        Triangulated panel whose every face normal is +Z.
    """
    if n_u < 1 or n_v < 1:
        raise ValueError(f"n_u and n_v must be >= 1, got {n_u}, {n_v}")
    if length_u <= 0.0 or width_v <= 0.0:
        raise ValueError(
            f"panel dimensions must be positive, got {length_u} x {width_v}")

    from surfmesh import convert_2d_face_to_3d, quad_faces_from_edges

    center = np.asarray(center, dtype=float)
    u_edges = np.linspace(-length_u / 2.0, length_u / 2.0, n_u + 1)
    v_edges = np.linspace(-width_v / 2.0, width_v / 2.0, n_v + 1)

    # surfmesh emits counter-clockwise quads in the (u, v) plane; lifting
    # them to Z = center_z keeps that winding, so the triangle normals come
    # out along +Z without any post-hoc flipping.
    quads = convert_2d_face_to_3d(quad_faces_from_edges(u_edges, v_edges),
                                  axis=2, offset=float(center[2]))
    quads[:, :, 0] += center[0]
    quads[:, :, 1] += center[1]

    data = np.zeros(2 * len(quads), dtype=mesh.Mesh.dtype)
    # Split each counter-clockwise quad (p0, p1, p2, p3) into the triangles
    # (p0, p1, p2) and (p0, p2, p3), preserving the winding.
    data['vectors'][0::2] = quads[:, [0, 1, 2], :]
    data['vectors'][1::2] = quads[:, [0, 2, 3], :]

    panel = mesh.Mesh(data)
    panel.update_normals()
    unit_normals = panel.get_unit_normals()
    assert np.allclose(unit_normals, [0.0, 0.0, 1.0], atol=1e-9), (
        'face normals must all point toward the plume source (+Z); the '
        'strike pipeline only loads faces whose normal opposes the plume')
    return panel


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--length-u', type=float, default=PANEL_LENGTH_U,
                        help='panel length along u / global X (m)')
    parser.add_argument('--width-v', type=float, default=PANEL_WIDTH_V,
                        help='panel width along v / global Y (m)')
    parser.add_argument('--n-u', type=int, default=DEFAULT_N_U,
                        help='quad divisions along u')
    parser.add_argument('--n-v', type=int, default=DEFAULT_N_V,
                        help='quad divisions along v')
    parser.add_argument('--out', type=str, default='iss_panel.stl')
    args = parser.parse_args()

    panel = build_panel_mesh(length_u=args.length_u, width_v=args.width_v,
                             n_u=args.n_u, n_v=args.n_v)
    out_path = Path(__file__).resolve().parent / args.out
    panel.save(str(out_path))
    print(f'saved {out_path} ({len(panel.vectors)} faces, '
          f'{args.length_u} m x {args.width_v} m, '
          f'{args.n_u} x {args.n_v} quads)')
