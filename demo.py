import numpy as np
import trimesh
# from MeshClass import MeshClass
# from smoothVF import smooth_vf
# from rdg_ADMM import rdg_ADMM
import matplotlib.pyplot as plt


# Define the filename and source point
filename = 'data/spot_rr.off'
x0 = [2273]  # source point / set

# # Load the mesh
mesh = trimesh.load(filename)
vertices = mesh.vertices
faces = mesh.faces

print(f"Type of vertices: {type(vertices)}")
print(f"Shape of vertices: {vertices.shape}")
print(f"Type of faces: {type(faces)}")
print(f"Shape of faces: {faces.shape}")

# Mm = MeshClass(vertices=vertices, faces=faces)

# print("MESH DONE")

# # Regularized - Dirichlet Energy
# alpha_hat0 = 0.05  # scale invariant, represents the weight of the regularizer

# u_D1 = rdg_ADMM(Mm, x0, reg='D', alpha_hat=alpha_hat0)
# u0 = rdg_ADMM(Mm, x0, reg='D', alpha_hat=0)  # No regularization
# u_D2 = rdg_ADMM(Mm, x0, reg='D', alpha_hat=3*alpha_hat0)  # Higher Regularization - Dirichlet Energy

# # Given directions
# given_vf_faces = [4736, 2703]
# given_vf_vals = np.array([[1.6256, -0.3518, -0.6234], [1.6952, 0.3193, 0.0335]])

# vf = np.zeros((Mm.nf, 3))
# vf[given_vf_faces, :] = given_vf_vals

# # Interpolate vf to mesh
# vf_int = smooth_vf(Mm, vf, 2)

# # Optionally, scale the interpolated line field with a geodesic Gaussian
# localize_vf = True
# if localize_vf:
#     vf_faces_v = Mm.faces[given_vf_faces, :].ravel()
#     dist_to_vf_faces = rdg_ADMM(Mm, vf_faces_v, reg='D', alpha_hat=0)
#     sigma2 = np.sum(Mm.ta) / 10**2
#     dist_vf_gaus = np.exp(-dist_to_vf_faces**2 / (2 * sigma2))
#     vf_int = Mm.interpulateVertices2Face(dist_vf_gaus)[0] * vf_int

# # Regularizers
# alpha_hat = 0.05
# beta_hat = 100

# u_vfa = rdg_ADMM(Mm, x0, reg='vfa', alpha_hat=alpha_hat, beta_hat=beta_hat, vf=vf_int)

# # Load camera settings
# cam = np.load('data/spot_rr_cam.npy', allow_pickle=True).item()['cam']

# u_all = np.concatenate([u0, u_D1, u_D2, u_vfa])
# umin = np.min(u_all)
# umax = np.max(u_all)
# nlines = 15

# # Visualize distances
# Mm.visualizeDistances(u0, x0, nlines, [umin, umax], cam)
# Mm.visualizeDistances(u_D1, x0, nlines, [umin, umax], cam)
# Mm.visualizeDistances(u_D2, x0, nlines, [umin, umax], cam)

# # Visualize vector field alignment
# Mm.visualizeDistances(u_vfa, x0, nlines, [umin, umax], cam)
# br = Mm.baryCentersCalc()
# plt.quiver(br[:, 0], br[:, 1], br[:, 2], vf[:, 0], vf[:, 1], vf[:, 2], length=2, color='k', linewidth=2)
# plt.quiver(br[:, 0], br[:, 1], br[:, 2], -vf[:, 0], -vf[:, 1], -vf[:, 2], length=2, color='k', linewidth=2)
# plt.show()
