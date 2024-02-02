from stl import mesh
from math import pi
plate = mesh.Mesh.from_file('flat_plate.stl')

plate.translate([-20, 0, 0])
theta = (pi/180) * 90
plate.rotate([1, 0, 0], theta = theta)
plate.translate([0, -7.0, 0])


plate.save('flate_plate_transformed.stl')
print(plate)