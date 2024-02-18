from stl import mesh
import math

mesh = mesh.Mesh.from_file('cluster_ATV216.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([-0.153,-0.1625,-0.18995]) # center on origin

mesh.translate([-0.153,0,0]) # line up edge of tc with edge of lm

mesh.save('cluster_ATV216_transformed.stl')