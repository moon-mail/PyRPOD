from stl import mesh
import numpy as np

# Load the STL file
your_mesh = mesh.Mesh.from_file('flat_plate_3.stl')

# Define the translation vector
translation_vector = np.array([-20, 0, 0])

# Translate the mesh
your_mesh.translate(translation_vector)

# Save the translated mesh to a new file
your_mesh.save('flat_plate_2.stl')
