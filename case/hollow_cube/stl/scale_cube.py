from stl import mesh
cube = mesh.Mesh.from_file('/home/nickpalumbo8/PyRPOD/data/stl/HollowCube.stl')
cube.points = 0.0025 * cube.points
cube.save('HollowCube_scale.stl')