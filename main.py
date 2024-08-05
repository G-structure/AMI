import potpourri3d as pp3d
import matplotlib.pyplot as plt
import numpy as np
import amigo

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

def normalize_unit_cube(X):
  # Normalize a mesh into a unit cube
  center = X.mean(0)
  X = X - center
  maxi = X.max(0)
  mini = X.min(0)
  X = X / (maxi - mini).max()
  return X

# load and normalize the vertex positions
V, F = load_off_file(filename)

#V = normalize_unit_cube(np.array(V))
#V += np.random.rand(V.shape\[0], V.shape\[1]) * 1e-4 # get rid of singularities
#F = np.array(F)

# Input: A triangle mesh $M$, seed $s$, stitch width $w$
# initialize the amigurumi object
m = amigo.Amigurumi(V, F)
seed = 0
w = 0.001

# Mesh to Graph :// Section4
## Geometry $\mathcal{S},X_{\mathcal{G}}$ :// Section4.1
s, xg = amigo.Geometry(m, seed, w)
print(s)
## Connectivity $\mathcal{R},\mathcal{C}$ :// Section4.2
r, c = amigo.Connectivity(m, seed, w, s, xg)

# Graph to pattern :// Section5
## Graph to program :// Section5.1
program = amigo.Program(m, seed, w, s, xg, r, c)
## Program to pattern $P(\mathcal{G})$ :// Section5.2
pattern = amigo.pattern(m, seed, w, s, xg, r, c, program)

# Output: Embedded crochet graph $\mathcal{G}=(\mathcal{S},\mathcal{R}\cup\mathcal{C}),X_{\mathcal{G}}$, crochet pattern $P(\mathcal{G})$
print(pattern)
