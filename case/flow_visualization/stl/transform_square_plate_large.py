from stl import mesh
import math

mesh = mesh.Mesh.from_file('square_plate_large_fine_defaced_series.STL')

# Center on origin
mesh.translate([0, -24, -24])
# mesh.translate([-1/2, -24, -24])

# Mate with YZ plane at x = 10 m
# mesh.translate([10.5, 0, 0])

mesh.save('square_plate_large_fine_defaced_series_transformed.stl')