from stl import mesh
mesh = mesh.Mesh.from_file('thruster_cluster_scaled.stl')

mesh.save('thruster_cluster_scaled_translated.stl')