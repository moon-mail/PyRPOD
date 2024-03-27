from stl import mesh
import math

mesh = mesh.Mesh.from_file('mold_funnel.stl')

mesh.translate([-0.217368/2, -1.6/2, -1.6/2])

mesh.save('mold_funnel_transformed.stl')