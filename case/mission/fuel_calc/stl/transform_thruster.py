from stl import mesh
import numpy as np
import math

# for example
cant_angle = 0

mesh = mesh.Mesh.from_file('thruster_ATV216.stl')

# center on origin
mesh.translate([-0.19/2,-0.095/2,-0.095/2])

# SweepConfig needs to manipulate the dcm and the thruster should be rotated before the following translation
# center point of thruster that intersects with cluster on origin
mesh.translate([0.095*np.cos(math.radians(cant_angle)),0,0.095*np.sin(math.radians(cant_angle))])

mesh.save('thruster_ATV216_transformed.stl')