include("MeshClass.jl")
include("smoothVF.jl")
include("RDG_ADMM.jl")

using .MeshClass
using .smoothVF
using .RDG_ADMM

begin
    using FileIO
    using GeometryTypes
    using GeometryBasics
    using MeshIO
    using GeometryBasics: coordinates, faces
end

# Clear variables, close all plots, and clear the console
# In Julia, you need to write equivalent code for your specific needs
# Example for clearing the console: Base.runtests(["tput", "reset"])
# Example for closing plots: Close relevant plot windows

# addpath("data/")

filename = "data/spot_rr.off"
x0 = 2273; # source point / set()

function to_int(x)
    return Int(GeometryBasics.value(x))
end

# load mesh()
mesh = load(filename)
# shape_vertices = [ [p[1], p[2], p[3]] for p in GeometryBasics.coordinates(mesh) ] |> hcat |> permutedims

function to_int(x)
    return Int(GeometryBasics.value(x))
end

faces_raw = GeometryBasics.faces(mesh)
# println("faces_raw: ", faces_raw)  # Debug: Check the structure of faces_raw

faces_int = [ [to_int(f[1]), to_int(f[2]), to_int(f[3])] for f in faces_raw ]
# println("faces_int: ", faces_int)  # Debug: Check the structure after conversion

face_matrix = hcat(faces_int...)   # The ... is crucial to splat the vectors
# println("face_matrix (hcat): ", face_matrix)  # Debug: Check the matrix after hcat

shape_faces = permutedims(face_matrix)
# println("face (permuted): ", faces)  # Debug: Check the final permuted face matrix

shape_vertices = reduce(hcat, [ [Float64(p[1]), Float64(p[2]), Float64(p[3])] for p in GeometryBasics.coordinates(mesh) ])
shape_vertices = reshape(collect(permutedims(shape_vertices)), (4764, 3))
println("Type of vertices: ", typeof(shape_vertices))
println("Shape of shape_vertices: ", size(shape_vertices))
println("Type of faces: ", typeof(shape_faces))
println("Shape of shape_faces: ", size(shape_faces))
Mm = meshClass(shape_vertices, shape_faces)

## Regularized - Dirichlet Energy
alpha_hat0 = 0.05; # scale invariant; represents the weight of the regularizer
                   # for Dirichlet regularizer - the size of the smoothing area()

u_D1 = rdg_ADMM[Mm, x0, "alpha_hat", alpha_hat0]

u0 = rdg_ADMM[Mm, x0, "alpha_hat", 0];               # No regularization
u_D2 = rdg_ADMM[Mm, x0, "alpha_hat", 3*alpha_hat0];  # Higher Regularization - Dirichlet Energy



## Regularized - Vector Field Alignment
# given directions:
given_vf_faces = [4736 2703]
given_vf_vals = [1.6256   -0.3518   -0.6234  1.6952    0.3193    0.0335]

vf = zeros(Mm.nf,3)
vf[given_vf_faces,:] = given_vf_vals

# interpolate vf to mesh()
vf_int = smooth_vf[Mm, vf, 2]

# Optionally; scale the interpolated line field with a geodesic Gaussian
localize_vf = 1
if localize_vf
    vf_faces_v = Mm.faces[given_vf_faces,:]; vf_faces_v = vf_faces_v[:]
    dist_to_vf_faces = rdg_ADMM[Mm, vf_faces_v, "alpha_hat", 0]
    sigma2 = sum(Mm.ta)/10^2; dist_vf_gaus = exp(-dist_to_vf_faces.^2/(2*sigma2))

    vf_int = Mm.interpulateVertices2Face[dist_vf_gaus].*vf_int
end


# regularizers
alpha_hat = 0.05
beta_hat = 100

u_vfa = rdg_ADMM[Mm, x0, "reg', 'vfa', 'alpha_hat', alpha_hat, 'beta_hat', beta_hat, 'vf", vf_int]



## Figures
cam = load("spot_rr_cam.mat"); cam = cam.cam
u_all = [u0[:] u_D1[:] u_D2[:] u_vfa[:]]'
umin = minimum(u_all)
umax = maximum(u_all)
nlines = 15

Mm.visualizeDistances[u0, x0, nlines, [umin, umax], cam]
Mm.visualizeDistances[u_D1, x0, nlines, [umin, umax], cam]
Mm.visualizeDistances[u_D2, x0, nlines, [umin, umax], cam]

Mm.visualizeDistances[u_vfa, x0, nlines, [umin, umax], cam]; hold on
br = Mm.baryCentersCalc
quiver3(br[:,1],br[:,2],br[:,3],...
    vf[:,1],vf[:,2],vf[:,3],2,"color','k','LineWidth',2,'ShowArrowHead','off")
quiver3(br[:,1],br[:,2],br[:,3],...
    -vf[:,1],-vf[:,2],-vf[:,3],2,"color','k','LineWidth',2,'ShowArrowHead','off")
