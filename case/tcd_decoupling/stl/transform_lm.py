from stl import mesh
import math

mesh = mesh.Mesh.from_file('logistics_module.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-12.5125,-5.317443,-5.317443])



mesh.save('logistics_module_transformed.stl')