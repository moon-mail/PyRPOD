from stl import mesh
import math

mesh = mesh.Mesh.from_file('square_plate_refinedx3.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-1/2, -12, -12])

mesh.save('square_plate_refinedx3_transformed.stl')