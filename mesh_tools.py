import argparse
import trimesh
import pymesh

def upsample_mesh(off_file_path, iterations=1):
    # Load the mesh from an OFF file
    mesh = trimesh.load(off_file_path)

    # Convert the Trimesh object to a PyMesh object
    vertices = mesh.vertices
    faces = mesh.faces
    pymesh_mesh = pymesh.form_mesh(vertices, faces)

    # Perform Loop subdivision
    for _ in range(iterations):
        pymesh_mesh, _ = pymesh.subdivide(pymesh_mesh, order=1, method="loop")

    # Convert back to a Trimesh object
    upsampled_mesh = trimesh.Trimesh(vertices=pymesh_mesh.vertices, faces=pymesh_mesh.faces)

    return upsampled_mesh

def main():
    parser = argparse.ArgumentParser(description="Upsample a triangular mesh from an OFF file.")
    parser.add_argument("input_file", type=str, help="Path to the input OFF file")
    parser.add_argument("output_file", type=str, help="Path to the output OFF file")
    parser.add_argument("-i", "--iterations", type=int, default=1, help="Number of upsampling iterations (default: 1)")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    iterations = args.iterations

    # Upsample the mesh with the specified number of iterations
    upsampled_mesh = upsample_mesh(input_file, iterations)

    # Export the upsampled mesh to an OFF file
    upsampled_mesh.export(output_file)

    print(f"Upsampled mesh saved to {output_file}")

if __name__ == "__main__":
    main()
