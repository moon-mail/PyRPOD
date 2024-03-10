from stl import mesh
import math

mesh = mesh.Mesh.from_file('flat_plate_defaced.STL')

mesh.points = 1000 * mesh.points

mesh.translate([-1/2, -70/2, -10/2])

mesh.save('flat_plate_defaced_transformed.STL')