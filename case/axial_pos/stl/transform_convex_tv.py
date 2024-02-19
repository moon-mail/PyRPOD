from stl import mesh

mesh = mesh.Mesh.from_file('convex_tv.stl')

mesh.points = 0.001 * mesh.points

mesh.translate([0,-56.533176/2,-56.533176/2])

mesh.save('convex_tv_transformed.stl')