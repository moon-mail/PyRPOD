from stl import mesh
from math import pi
plate = mesh.Mesh.from_file('square_plate_refined.stl')

plate.translate([-2.5, 0, -2.5])
theta = (pi/180) * 90
plate.rotate([1, 0, 0], theta = theta)
plate.points = 5 * plate.points


plate.save('square_plate_transformed.stl')
print(plate)