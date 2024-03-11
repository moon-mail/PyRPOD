from stl import mesh
import math

mesh = mesh.Mesh.from_file('square_plate_low_res_defaced.STL')

mesh.points = 1000 * mesh.points

mesh.translate([-1/2, -12, -12])

mesh.save('square_plate_low_res_defaced_transformed.STL')