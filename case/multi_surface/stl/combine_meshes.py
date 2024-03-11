from stl import Mesh

out_fh = open('multi_plate.stl', 'wb')
for i in range(2):
    iter = str(i + 1)
    mesh = Mesh.from_file('flat_plate_' + iter + '.stl')
    mesh.save('mesh_%d.stl' % i, out_fh)
out_fh.close()