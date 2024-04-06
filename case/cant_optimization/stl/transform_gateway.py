from stl import mesh
import math

mesh = mesh.Mesh.from_file('gateway_med_defaced.stl')

mesh.translate([-1.95672025, -23.00022212, -58.41826258])

mesh.save('gateway_med_defaced_transformed.stl')