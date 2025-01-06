from stl import mesh
import math

mesh = mesh.Mesh.from_file('cluster.stl')

mesh.translate([-0.306/2,-0.325/2,-0.3799/2]) # center on origin

mesh.translate([-0.306/2,0,0]) # line up edge of tc with edge of lm

mesh.save('cluster_transformed.stl')