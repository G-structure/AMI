import numpy as np
import potpourri3d as pp3d
import scipy.sparse as sparse
import scipy.sparse.linalg as spla

def get_boundary_vertices(cut_path, V):
    """
    Get the boundary vertices along the cut path.
    """
    # Assume cut_path is a list of vertex indices
    return np.unique(cut_path)

def gradient_operator(V, F):
    num_vertices = V.shape[0]
    num_faces = F.shape[0]

    # Compute face normals and areas
    v1 = V[F[:, 1]] - V[F[:, 0]]
    v2 = V[F[:, 2]] - V[F[:, 0]]
    face_normals = np.cross(v1, v2)
    face_areas = np.linalg.norm(face_normals, axis=1) / 2
    face_normals /= np.linalg.norm(face_normals, axis=1)[:, np.newaxis]

    # Compute per-face gradients
    G = np.zeros((num_faces, 3, 3))
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        eij = V[F[:, j]] - V[F[:, i]]
        G[:, i] = np.cross(face_normals, eij) / (2 * face_areas[:, np.newaxis])

    # Construct sparse gradient operator
    I = np.repeat(np.arange(num_faces), 3)
    J = F.flatten()
    V = G.reshape(-1, 3)

    return sparse.coo_matrix((V.flatten(), (I.repeat(3), np.column_stack([J]*3).flatten())), shape=(num_faces, num_vertices))

def compute_Jf(f, G, V, F):
    num_vertices = V.shape[0]
    num_faces = F.shape[0]

    # Compute the gradient of f
    grad_f = G @ f

    # Reshape grad_f to match the number of faces
    grad_f = grad_f.reshape(num_faces, 3)

    # Compute face normals
    v1 = V[F[:, 1]] - V[F[:, 0]]
    v2 = V[F[:, 2]] - V[F[:, 0]]
    face_normals = np.cross(v1, v2)
    face_normals /= np.linalg.norm(face_normals, axis=1)[:, np.newaxis]

    # Rotate gradient by π/2 in the tangent plane
    Jgrad_f = np.cross(face_normals, grad_f)

    # Average rotated gradients to vertices
    Jf = np.zeros((num_vertices, 3))
    np.add.at(Jf, F[:, 0], Jgrad_f)
    np.add.at(Jf, F[:, 1], Jgrad_f)
    np.add.at(Jf, F[:, 2], Jgrad_f)

    # Normalize, avoiding division by zero
    norms = np.linalg.norm(Jf, axis=1)
    Jf[norms > 0] /= norms[norms > 0, np.newaxis]

    return Jf.flatten()

def Geometry(amigo, seed_idx, w):
    # Compute geodesic distance from seed (f function)
    f = amigo.solver.compute_distance(seed_idx)

    # Compute g (orthogonal to f)
    # We need to solve the optimization problem to find g

    # First, create the cut along the geodesic from seed to fM
    fM_idx = np.argmax(f)
    # Compute the direction from seed to fM
    direction = amigo.V[fM_idx] - amigo.V[seed_idx]
    direction_normalized = direction / np.linalg.norm(direction)

    # Trace the geodesic from seed to fM
    trace_pts = amigo.tracer.trace_geodesic_from_vertex(seed_idx, direction)

    # Create the cut along this geodesic path
    cut_path = np.array([np.argmin(np.linalg.norm(amigo.V - pt, axis=1)) for pt in trace_pts])

    # Create Laplacian matrix
    L = pp3d.cotan_laplacian(amigo.V, amigo.F)

    # Create the gradient operator
    G = gradient_operator(amigo.V, amigo.F)

    # Compute J∇f
    Jf = compute_Jf(f, G, amigo.V, amigo.F)

    # Set up the linear system to solve for g
    A = L.tocsr()
    b = G.T @ Jf.reshape(-1, 1)

    # Add boundary conditions
    boundary_vertices = get_boundary_vertices(cut_path, amigo.V)
    for v in boundary_vertices:
        A[v, :] = 0
        A[v, v] = 1
        b[v] = 0

    # Solve the system
    g = spla.spsolve(A, b)

    # Sample f and g
    f_samples = np.arange(0, f.max(), w)
    g_samples = np.arange(g.min(), g.max(), w)

    S = []
    X_G = []
    for f_val in f_samples:
        for g_val in g_samples:
            idx = np.argmin((f - f_val)**2 + (g - g_val)**2)
            S.append(idx)
            X_G.append(amigo.V[idx])

    return np.array(S), np.array(X_G)

class Amigurumi:
    def __init__(self, V, F):
        self.V = np.array(V)
        self.F = np.array(F)
        self.solver = pp3d.MeshHeatMethodDistanceSolver(self.V, self.F)
        self.vector_solver = pp3d.MeshVectorHeatSolver(self.V, self.F)
        self.tracer = pp3d.GeodesicTracer(self.V, self.F)

    # def Geometry(self, seed, w):
    #     # Compute geodesic distance from seed
    #     f = self.solver.compute_distance(seed)

    #     # Compute g (orthogonal to f)
    #     basisX, basisY, _ = self.vector_solver.get_tangent_frames()
    #     g = self.vector_solver.extend_scalar([seed], [0.0])

    #     # Sample f and g
    #     f_samples = np.arange(0, f.max(), w)
    #     # Handle the case where g range is too small
    #     g_min, g_max = g.min(), g.max()
    #     if g_min == g_max:
    #         g_samples = np.array([g_min])
    #     else:
    #         g_range = g_max - g_min
    #         num_g_samples = max(int(g_range / w), 1)  # Ensure at least one sample
    #         g_samples = np.linspace(g_min, g_max, num_g_samples)

    #     for value in g_samples:
    #         print(f"hi {value}")

    #     S = []
    #     X_G = []
    #     for f_val in f_samples:
    #         for g_val in g_samples:
    #             idx = np.argmin((f - f_val)**2 + (g - g_val)**2)
    #             S.append(idx)
    #             X_G.append(self.V[idx])

    #     return np.array(S), np.array(X_G)

def Connectivity(amigo, seed, w, S, X_G):
    R = []  # Row edges
    C = []  # Column edges

    for i in range(len(S) - 1):
        if f[S[i]] == f[S[i+1]]:
            R.append((S[i], S[i+1]))
        else:
            C.append((S[i], S[i+1]))

    return np.array(R), np.array(C)

def Program(amigo, seed, w, S, X_G, R, C):
    program = []
    current_row = 0
    stitches_in_row = 1

    for i, s in enumerate(S):
        if i > 0 and f[S[i]] != f[S[i-1]]:
            program.append(f"Row {current_row}: {stitches_in_row} sc")
            current_row += 1
            stitches_in_row = 0

        if (s, S[i+1]) in R:
            stitches_in_row += 1
        elif (s, S[i+1]) in C:
            if stitches_in_row > 1:
                program.append("inc")
            else:
                program.append("dec")

    program.append(f"Row {current_row}: {stitches_in_row} sc")
    return program

def pattern(amigo, seed, w, S, X_G, R, C, program):
    # Simplify the program into a human-readable pattern
    pattern = []
    current_instruction = ""
    count = 0

    for instruction in program:
        if instruction == current_instruction:
            count += 1
        else:
            if count > 0:
                pattern.append(f"{count}x {current_instruction}")
            current_instruction = instruction
            count = 1

    if count > 0:
        pattern.append(f"{count}x {current_instruction}")

    return pattern
