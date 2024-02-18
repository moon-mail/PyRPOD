from stl import mesh
import math

mesh = mesh.Mesh.from_file('mold_funnel.stl')

mesh.translate([0, 0, -54.342])
# Additional translations used for mold_funnel_centerline, 1000 is the length of the centerline and 39.728 is exit diameter
# plumeMesh.translate([-39.728/2, -39.728/2, -1000])    # nozzle throat (x, y, z) = (0, 0, 0)
mesh.rotate([1,0,0], math.radians(180))

mesh.points = 0.05 * mesh.points

mesh.save('mold_funnel_transformed.stl')