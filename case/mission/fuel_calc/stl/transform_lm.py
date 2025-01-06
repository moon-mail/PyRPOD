from stl import mesh
import math

mesh = mesh.Mesh.from_file('logistics_module.stl')

mesh.points = 1000 * mesh.points

# center radially on x-axis and align face of LM's docking port to plane at x = 0
mesh.translate([-12.105,-9.659079/2,-9.659079/2])

# align edge of LM's outer diameter on the docking side to plane at x = 0
mesh.translate([0.46,0,0])

mesh.save('logistics_module_transformed.stl')