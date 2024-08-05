import redis
import time
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize

r = redis.Redis(host='localhost', port=6379, db=0)

queue_name = 'matlab_data'

def generate_crochet_pattern(vertices, faces, distances):
    instructions = []
    seed_vertex = int(np.argmin(distances))  # Convert NumPy int64 to native int
    max_vertex = int(np.argmax(distances))
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

    pattern = {
        'seed_vertex': seed_vertex,
        'row_edges': row_edges,
        'col_edges': col_edges,
        'stitches': instructions
    }
    return pattern

def find_row_order(vertices, faces, distances, w):
    max_distance = np.max(distances)
    isoline_values = np.arange(0, max_distance, w)
    row_order = []

    for isoline_value in isoline_values:
        close_vertices = vertices[np.abs(distances - isoline_value) < w / 2]
        row = [
            {
                'vert_index': int(i),  # Convert NumPy int64 to native int
                'raw_distance': float(distances[i]),  # Convert NumPy float to native float
                'reg_distance': isoline_value
            }
            for i, vertex in enumerate(vertices)
            if np.abs(distances[i] - isoline_value) < w / 2
        ]
        row_order.append(row)

    return row_order


def compute_gradient(vertices, values):
    gradients = np.zeros_like(vertices)
    for i, v in enumerate(vertices):
        neighbors = np.where(np.linalg.norm(vertices - v, axis=1) < 1e-5)[0]
        if len(neighbors) > 1:
            deltas = vertices[neighbors] - v
            value_deltas = values[neighbors] - values[i]
            gradients[i] = np.mean(value_deltas[:, np.newaxis] * deltas, axis=0)
    return gradients

def rotated_gradient(gradients):
    J = np.array([[0, -1], [1, 0]])  # Rotation matrix for 2D case
    return np.dot(gradients[:, :2], J.T)

def compute_gradient_g(g, vertices):
    gradient_g = np.zeros_like(vertices)
    for i, v in enumerate(vertices):
        neighbors = np.where(np.linalg.norm(vertices - v, axis=1) < 1e-5)[0]
        if len(neighbors) > 1:
            deltas = vertices[neighbors] - v
            g_deltas = g[neighbors] - g[i]
            gradient_g[i] = np.mean(g_deltas[:, np.newaxis] * deltas, axis=0)
    return gradient_g

def objective(g, vertices, gradients, rotated_gradients):
    gradient_g = compute_gradient_g(g, vertices)
    diff = np.einsum('ij,ij->i', gradients, gradient_g) - 1
    return np.sum(diff**2)

def find_column_order(vertices, f_values, boundary_indices):
    gradients = compute_gradient(vertices, f_values)
    print(gradients)
    rotated_gradients = rotated_gradient(gradients)

    g_initial = np.zeros(len(vertices))
    g_initial[boundary_indices] = 0

    result = minimize(objective, g_initial, args=(vertices, gradients, rotated_gradients), method='L-BFGS-B')

    if not result.success:
        print("Optimization failed:", result.message)

    return result.x

def compute_geodesic_path(vertices, faces, start_idx, end_idx):
    n_vertices = len(vertices)
    adj_matrix = np.zeros((n_vertices, n_vertices))

    for face in faces:
        for i in range(3):
            for j in range(i + 1, 3):
                v0, v1 = face[i], face[j]
                distance = np.linalg.norm(vertices[v0] - vertices[v1])
                adj_matrix[v0, v1] = distance
                adj_matrix[v1, v0] = distance

    graph = csr_matrix(adj_matrix)
    dist_matrix, predecessors = shortest_path(csgraph=graph, directed=False, indices=start_idx, return_predecessors=True)

    path = []
    current_vertex = end_idx
    while current_vertex != start_idx:
        path.append(current_vertex)
        current_vertex = predecessors[current_vertex]
    path.append(start_idx)
    path.reverse()

    return path

def cut_mesh(vertices, faces, geodesic_path):
    new_vertices = vertices.tolist()
    new_faces = []
    vertex_map = {}

    for v in geodesic_path:
        new_idx = len(new_vertices)
        vertex_map[v] = new_idx
        new_vertices.append(vertices[v])

    for face in faces:
        new_face = []
        for v in face:
            if v in geodesic_path:
                new_face.append(vertex_map[v])
            else:
                new_face.append(v)
        new_faces.append(new_face)

    new_vertices = np.array(new_vertices)
    new_faces = np.array(new_faces)

    return new_vertices, new_faces

while True:
    item = r.lpop(queue_name)
    if item:
        data = json.loads(item.decode('utf-8'))
        vertices = np.array(data['vertices'])
        faces = np.array(data['faces'])
        distances = np.array(data['u_vfaOut'])

        row_order = find_row_order(vertices, faces, distances, 0.03)

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.azim = 90  # Rotate the plot to match the desired view
        # ax.elev = -90
        # ax.roll = 45

        for close_vertices in row_order:
            vertices_to_plot = np.array([vertices[item['vert_index']] for item in close_vertices])
            # ax.scatter(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2], s=1)

        plt.figure()
        column_orders = []
        for close_vertices in row_order:
            if len(close_vertices) > 1:
                f_values = np.array([item['reg_distance'] for item in close_vertices])
                row_vertices = np.array([vertices[item['vert_index']] for item in close_vertices])
                boundary_indices = [0]  # Example boundary index, replace with actual boundary condition
                g_values = find_column_order(row_vertices, f_values, boundary_indices)
                column_orders.append(g_values)
                plt.scatter(g_values, f_values, c='blue', s=10)

        plt.xlabel('g values')
        plt.ylabel('f values')
        plt.title('Row Order by g values')
        plt.grid(True)
        plt.show()

        instructions = generate_crochet_pattern(vertices, faces, distances)
        with open('crochet_instructions.json', 'w') as f:
            json.dump(instructions, f, indent=2)

        # plt.show()
        print(f"Processed data item from the queue.")
    else:
        print("Queue is empty, waiting for new items...")
        time.sleep(2)
