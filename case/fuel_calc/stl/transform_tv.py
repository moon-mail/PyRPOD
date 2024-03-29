from stl import mesh
import math

mesh = mesh.Mesh.from_file('convex_tv.stl')

# convex_tv
mesh.translate([0,-56.5331761/2,-56.5331761/2])

mesh.save('convex_tv_transformed.stl')