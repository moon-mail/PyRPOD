from stl import mesh
mesh = mesh.Mesh.from_file('thruster_cluster.stl')
mesh.points = 20 * mesh.points
mesh.save('thruster_cluster_scaled.stl')