from stl import mesh
import math

mesh = mesh.Mesh.from_file('convex_tv.stl')

mesh.translate([0, -27.5, -27.5])

mesh.save('convex_tv_transformed.stl')