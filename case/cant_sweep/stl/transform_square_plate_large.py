from stl import mesh
import math

mesh = mesh.Mesh.from_file('square_plate_large.stl')

mesh.translate([-1/2, -24, -24])

mesh.save('square_plate_large_transformed.stl')