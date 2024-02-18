import numpy as np
from stl import mesh
from tqdm import tqdm

def subdivide_triangle(p1, p2, p3):
    # Calculate midpoints of each side
    m1 = (p1 + p2) / 2
    m2 = (p2 + p3) / 2
    m3 = (p1 + p3) / 2

    # Return the four new triangles
    return [
        (p1, m1, m3),
        (p2, m2, m1),
        (p3, m3, m2),
        (m1, m2, m3),
    ]

def increase_mesh_fidelity(input_file, output_file, refinement_levels=2):
    # Load the mesh
    original_mesh = mesh.Mesh.from_file(input_file)
    new_triangles = []

    # Wrap original_mesh with tqdm to display a progress bar
    for i in tqdm(range(len(original_mesh)), desc="Refining Mesh"):
        v1, v2, v3 = original_mesh.vectors[i]
        # Start with initial triangles
        triangles_to_refine = [(v1, v2, v3)]

        # Refine each triangle by the specified number of refinement levels
        for _ in range(refinement_levels):
            refined_triangles = []
            for triangle in triangles_to_refine:
                # Subdivide each triangle and add the new ones to the list to refine in the next iteration
                refined_triangles.extend(subdivide_triangle(*triangle))
            triangles_to_refine = refined_triangles

        # After all refinement levels, add the final set of triangles to the new_triangles list
        new_triangles.extend(triangles_to_refine)

    # Create a new mesh object
    new_mesh = mesh.Mesh(np.zeros(len(new_triangles), dtype=mesh.Mesh.dtype))

    # Fill the mesh with new triangles
    for i, triangle in enumerate(new_triangles):
        new_mesh.vectors[i] = np.array(triangle)

    # Save the new mesh
    new_mesh.save(output_file)
    print(f"New mesh with increased fidelity saved to '{output_file}'")

# Run
input_stl = 'hollow_cube_transformed_refined.stl'
output_stl = 'hollow_cube_transformed_refinedx2.stl'
increase_mesh_fidelity(input_stl, output_stl)