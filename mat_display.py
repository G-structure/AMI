import redis
import time
import json
import numpy as np  # Import numpy for array manipulations
import matplotlib.pyplot as plt  # Import matplotlib for plotting
from mpl_toolkits.mplot3d import Axes3D  # Import module for 3D plotting

r = redis.Redis(host='localhost', port=6379, db=0)

queue_name = 'matlab_data'
def generate_crochet_pattern(vertices, faces, distances):
    """
    Generate crochet pattern instructions from 3D mesh data.
    """
    instructions = []
    seed_vertex = int(np.argmin(distances))  # Convert NumPy int64 to native int
    row_edges = []
    col_edges = []
    S = [[]]



    # Step 1: Select seed and calculate geodesic distances (already done)
    # Step 2: Iterate vertices and create rows of stitches
    for i, vertex in enumerate(vertices):
        stitch = {
            'vertex_index': int(i),  # Convert NumPy int64 to native int
            'position': vertex.tolist(),
            'distance': float(distances[i])  # Convert NumPy float to native float
        }
        if i > 0:
            row_edges.append((int(i - 1), int(i)))  # Convert NumPy int64 to native int
        instructions.append(stitch)

    # Step 3: Form columns (inter-row connections)
    for start, end in zip(row_edges[:-1], row_edges[1:]):
        col_edges.append((int(start[1]), int(end[0])))  # Convert NumPy int64 to native int

    # Combine into instructions
    pattern = {
        'seed_vertex': seed_vertex,
        'row_edges': row_edges,
        'col_edges': col_edges,
        'stitches': instructions
    }
    return pattern

while True:
    item = r.lpop(queue_name)
    if item:
        data = json.loads(item.decode('utf-8'))
        vertices = np.array(data['vertices'])
        faces = np.array(data['faces'])
        distances = np.array(data['u_vfaOut'])

        # Visualize u_vfaOut
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c=distances, cmap='viridis')
        fig.colorbar(scatter, ax=ax, label='Geodesic Distance')
        # plt.show()

        # Generate crochet pattern from vertices and distances
        instructions = generate_crochet_pattern(vertices, faces, distances)
        # print(f"Crochet Instructions: {json.dumps(instructions, indent=2)}")
        with open('crochet_instructions.json', 'w') as f:
            json.dump(instructions, f, indent=2)

        print(f"Processed data item from the queue.")
    else:
        print("Queue is empty, waiting for new items...")
        time.sleep(2)  # Wait for 2 seconds before checking again
