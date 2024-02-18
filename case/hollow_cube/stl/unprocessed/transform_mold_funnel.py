from stl import mesh
import math

mesh = mesh.Mesh.from_file('mold_funnel.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-0.054342*4/2, -0.2*4, -0.2*4])

mesh.save('mold_funnel_transformed.stl')