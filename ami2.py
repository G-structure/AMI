import potpourri3d as pp3d
import matplotlib.pyplot as plt
import numpy as np

# filename = 'data/bent_pipe_closed_lr.off'
filename = 'data/spot_rr.off'

def load_off_file(filename):
    with open(filename, 'r') as file:
        if 'OFF' != file.readline().strip():
            raise ValueError('Not a valid OFF header')

        n_verts, n_faces, _ = map(int, file.readline().strip().split(' '))

        vertices = []
        for _ in range(n_verts):
            vertices.extend(map(float, file.readline().strip().split(' ')))

        vertices = [vertices[i:i + 3] for i in range(0, len(vertices), 3)]

        faces = []
        for _ in range(n_faces):
            face = list(map(int, file.readline().strip().split(' ')))
            faces.append(face[1:])

    return vertices, faces

def find_row_order(vertices, faces, distances, w):
    max_distance = np.max(distances)
    min_distance = np.min(distances)
    isoline_values = np.arange(min_distance, max_distance, w)
    row_order = []

    row_order.append([{
        'vert_index': int(np.argmax(distances)),
        'raw_distance': float(min_distance),
        'reg_distance': min_distance
    }])

    for isoline_value in isoline_values[1:-1]:
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

    row_order.append([{
        'vert_index': int(np.argmax(distances)),
        'raw_distance': float(max_distance),
        'reg_distance': max_distance
    }])

    return row_order

V, F = load_off_file(filename)
V = np.array(V)
F = np.array(F)
solver = pp3d.MeshHeatMethodDistanceSolver(V, F)
dist = solver.compute_distance(0)

row_order = find_row_order(V, F, dist, 0.03)

path_solver = pp3d.EdgeFlipGeodesicSolver(V,F)

col_order = []
# Append first row
col_order.append([item['vert_index'] for item in row_order[0]])
for row in row_order[1:-1]:
    v_indexs = [item['vert_index'] for item in row]
    if len(v_indexs) > 1:
        out = path_solver.find_geodesic_path_poly(v_indexs)
        matching_indices = [i for i in v_indexs if V[i].tolist() in out]
        # print(matching_indices)
        col_order.append(matching_indices)
# Append last row
col_order.append([item['vert_index'] for item in row_order[-1]])
S = col_order


# Create a 3D scatter plot of the vertices, colored by the geodesic distances
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(np.array(V)[:, 0], np.array(V)[:, 1], np.array(V)[:, 2], c=dist, cmap='viridis', s=10)
fig.colorbar(scatter, ax=ax, label='Geodesic Distance')
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')
ax.set_title('Vertices Colored by Geodesic Distances')

fig2, ax2 = plt.subplots()
for i, row in enumerate(S):
    for j, point in enumerate(row):
        ax2.scatter(j, i, color='r', s=10)
ax2.set_xlabel('g')
ax2.set_ylabel('f')
ax2.set_title('The f , g parameterization, and the sampled grid S')
plt.show()
