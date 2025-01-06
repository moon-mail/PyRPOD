from stl import mesh
import math

mesh = mesh.Mesh.from_file('square_plate_large_med.stl')

# Center on origin
mesh.translate([-1/2, -24, -24])

# Mate with YZ plane at x = 5 m
mesh.translate([5.5, 0, 0])

mesh.save('square_plate_large_med_transformed.stl')