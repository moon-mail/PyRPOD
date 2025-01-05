from stl import mesh
import math

mesh = mesh.Mesh.from_file('cylinder.STL')

mesh.points = 0.001 * mesh.points

mesh.translate([-7, -2, -2])

mesh.save('cylinder_transformed.STL')