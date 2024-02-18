from stl import mesh
import math

mesh = mesh.Mesh.from_file('hollow_cube.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-75/2, -106.066017/2, -106.066017/2])

mesh.save('hollow_cube_transformed.stl')