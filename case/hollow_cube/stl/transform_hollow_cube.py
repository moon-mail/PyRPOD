from stl import mesh
import math

mesh = mesh.Mesh.from_file('hollow_cube_defaced.STL')

mesh.points = 1000 * mesh.points

mesh.translate([-75/2, -106.066017/2, -106.066017/2])

mesh.save('hollow_cube_defaced_transformed.STL')