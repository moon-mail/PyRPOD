import argparse
import json
from stl import mesh
import numpy as np

def transform_mesh(input_file, output_file, scale, translate):
    # Load the mesh from the file
    model = mesh.Mesh.from_file(input_file)

    # Scale the mesh
    model.points *= scale

    # Translate the mesh
    model.translate(translate)

    # Save the transformed mesh to the output file
    model.save(output_file)
    print(f"Mesh saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform an STL mesh with scaling and translation using a configuration file.")

    # Configuration file
    parser.add_argument("config_file", type=str, help="Path to the configuration JSON file.")

    args = parser.parse_args()

    # Load configuration from the file
    with open(args.config_file, 'r') as config_file:
        config = json.load(config_file)

    input_file = config.get("input_file")
    output_file = config.get("output_file")
    scale = config.get("scale", 1.0)
    translate = config.get("translate", [0.0, 0.0, 0.0])

    # Transform the mesh using the configuration
    transform_mesh(
        input_file=input_file,
        output_file=output_file,
        scale=scale,
        translate=translate
    )
