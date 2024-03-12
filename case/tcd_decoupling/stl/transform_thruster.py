from stl import mesh
import numpy as np
import math

cant_angle = 0

mesh = mesh.Mesh.from_file('thruster_ATV216.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-0.095,-0.0475,-0.0475]) # center on origin

mesh.translate([0.095*np.cos(math.radians(cant_angle)),0,0.095*np.sin(math.radians(cant_angle))]) # center point of thruster that intersects with cluster on origin

mesh.save('thruster_ATV216_transformed.stl')