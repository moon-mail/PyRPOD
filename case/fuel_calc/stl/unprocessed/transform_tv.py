from stl import mesh
import math

mesh = mesh.Mesh.from_file('convex_tv.stl')

# center on origin
# mesh.translate([-64,-64,-64]) # KoS

# convex_tv
mesh.translate([16,-56.533176/2,-56.533176/2])

mesh.save('convex_tv_transformed.stl')