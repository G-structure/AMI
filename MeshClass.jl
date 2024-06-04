module MeshClass

using SparseArrays
using LinearAlgebra

include("CotLaplacian.jl")
using .CotLaplacian

export meshClass, normalize_vf, normv, get_camera, set_camera

# Define the meshClass struct
mutable struct meshClass
    vertices::Matrix{Float64}
    faces::Matrix{Int64}
    edges::Matrix{Int64}

    nv::Int64
    nf::Int64
    ne::Int64
    nie::Int64

    va::Vector{Float64}
    ta::Vector{Float64}
    ea::Vector{Float64}

    Nf::Matrix{Float64}
    Nv::Matrix{Float64}

    E1::Matrix{Float64}
    E2::Matrix{Float64}
    E3::Matrix{Float64}
    E::Matrix{Int64}
    R::SparseMatrixCSC{Float64,Int64}

    F1::Matrix{Float64}
    F2::Matrix{Float64}
    EB::SparseMatrixCSC{Float64,Int64}
    EBI::SparseMatrixCSC{Float64,Int64}

    e2t::Matrix{Int64}
    t2e::Matrix{Int64}
    v2e::SparseMatrixCSC{Int64,Int64}
    e2t1::SparseMatrixCSC{Float64,Int64}
    e2t2::SparseMatrixCSC{Float64,Int64}

    ie::Vector{Int64}
    inner_edges::Vector{Int64}
    bv::Vector{Int64}
    bf::Vector{Int64}

    Lap::SparseMatrixCSC{Float64,Int64}
    Ww::SparseMatrixCSC{Float64,Int64}
    Aa::Diagonal{Float64,Vector{Float64}}

    G::SparseMatrixCSC{Float64,Int64}
    D::SparseMatrixCSC{Float64,Int64}
end

function meshClass(vertices, faces)
    obj = meshClass(
        vertices, faces, zeros(Int64, 1, 1),
        size(vertices, 1), size(faces, 1), 0, 0,
        zeros(0), zeros(0), zeros(0),
        zeros(0, 3), zeros(0, 3),
        zeros(0, 3), zeros(0, 3), zeros(0, 3),
        zeros(0, 2), sparse(zeros(1, 1)),
        zeros(0, 3), zeros(0, 3), sparse(zeros(1, 1)), sparse(zeros(1, 1)),
        zeros(0, 4), zeros(0, 3), sparse(zeros(1, 1)),
        sparse(zeros(1, 1)), sparse(zeros(1, 1)),
        zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0),
        sparse(zeros(1, 1)), sparse(zeros(1, 1)), Diagonal(zeros(0)),
        sparse(zeros(1, 1)), sparse(zeros(1, 1))
    )
    compute_all!(obj)
    return obj
end

function compute_all!(obj::meshClass)
    obj.nv = size(obj.vertices, 1)
    obj.nf = size(obj.faces, 1)
    Nf = reduce(vcat, [cross(a, b)' for (a, b) in zip(eachrow(obj.vertices[obj.faces[:, 1], :] .- obj.vertices[obj.faces[:, 2], :]),
                                                      eachrow(obj.vertices[obj.faces[:, 1], :] .- obj.vertices[obj.faces[:, 3], :]))])
    obj.Nf = normalize_vf(Nf)
    obj.ta = vec(sqrt.(sum(Nf .^ 2, dims=2)) ./ 2)
    obj.Nv = vertex_normals(obj)
    println("Type of nv: ", typeof(obj.nv))
    println("Value of nv: ", obj.nv)

    println("Type of nf: ", typeof(obj.nf))
    println("Value of nf: ", obj.nf)

    println("Type of Nf: ", typeof(obj.Nf))
    println("Shape of Nf: ", size(obj.Nf))
    println("Value of Nf: ", obj.Nf[2,2])

    println("Type of ta: ", typeof(obj.ta))
    println("Shape of ta: ", size(obj.ta))
    println("Value of ta: ", obj.ta[2])

    println("Type of Nv: ", typeof(obj.Nv))
    println("Shape of Nv: ", size(obj.Nv))
    println("Value of Nv: ", obj.Nv[2,2])

    obj.va = calculatefvConnectivity(obj)' * obj.ta ./ 3
    obj.E1 = obj.vertices[obj.faces[:, 2], :] .- obj.vertices[obj.faces[:, 3], :]
    obj.E2 = obj.vertices[obj.faces[:, 1], :] .- obj.vertices[obj.faces[:, 3], :]
    obj.E3 = obj.vertices[obj.faces[:, 1], :] .- obj.vertices[obj.faces[:, 2], :]
    println("Type of Nv: ", typeof(obj.Nv))
    println("Shape of Nv: ", size(obj.Nv))
    println("Value of Nv: ", obj.Nv[2,2])
    obj.R = rot(obj)
    EE = sort([obj.faces[:, 1] obj.faces[:, 2]; obj.faces[:, 2] obj.faces[:, 3]; obj.faces[:, 3] obj.faces[:, 1]], dims=2)
    obj.E = unique(EE, dims=1)
    obj.edges, obj.e2t, obj.t2e, obj.e2t1, obj.e2t2, obj.v2e, obj.ie, obj.ne, obj.nie, obj.inner_edges, obj.bv, obj.bf = nc_data(obj)
    obj.ea = edge_areas(obj)
    edge_basis!(obj)
    compute_LB!(obj)
    GG!(obj)
    DD!(obj)
end

function calculatefvConnectivity(obj::meshClass)
    return sparse(vec(repeat(1:obj.nf, 1, 3)), vec(obj.faces), vec(ones(size(obj.faces))), obj.nf, obj.nv)
end

# function calculatefvConnectivity(obj::meshClass)
#     return sparse(repeat(1:obj.nf, 1, 3), vec(obj.faces), ones(size(obj.faces)), obj.nf, obj.nv)
# end

function edge_areas(obj::meshClass)
    T = convert(Matrix{Int64}, obj.faces)
    I = vec(hcat(T[:, 2], T[:, 3], T[:, 1]))
    J = vec(hcat(T[:, 3], T[:, 1], T[:, 2]))
    S = repeat(obj.ta / 3, inner=3)
    In = vcat(I, J)
    Jn = vcat(J, I)
    Sn = vcat(S, S)
    W = sparse(In, Jn, Sn, obj.nv, obj.nv)
    ea = zeros(length(obj.edges))
    for i in 1:length(ea)
        ea[i] = W[obj.edges[i, 1], obj.edges[i, 2]]
        ea[i] += W[obj.edges[i, 2], obj.edges[i, 1]]
    end
    return ea
end

function edge_basis!(mesh::meshClass)
    mesh.F1 = normalize_vf(mesh.E1)
    mesh.F2 = reshape(mesh.R * vec(mesh.F1), :, 3)

    I = repeat(1:mesh.nf, 1, 3)
    I = I[:]
    J = (1:mesh.nf)
    J = vcat(J, J + mesh.nf, J + 2 * mesh.nf)

    B1 = sparse(I, J, vec(mesh.F1), mesh.nf, 3 * mesh.nf)
    B2 = sparse(I, J, vec(mesh.F2), mesh.nf, 3 * mesh.nf)

    mesh.EB = hcat(B1, B2)'
    mesh.EBI = mesh.EB'
end

function rot(mesh::meshClass)
    sf = mesh.nf
    n = mesh.Nf

    II = repeat(1:sf, 1, 2)
    II = [II II .+ sf II .+ 2 * sf]'
    II = II[:]
    JJ1 = (1:sf)
    JJ2 = JJ1 .+ sf
    JJ3 = JJ2 .+ sf
    JJ = [JJ2 JJ3 JJ1 JJ3 JJ1 JJ2]'
    JJ = JJ[:]
    SS = [-n[:, 3] n[:, 2] n[:, 3] -n[:, 1] -n[:, 2] n[:, 1]]'
    SS = SS[:]

    return sparse(II, JJ, SS, 3 * sf, 3 * sf)
end

function vertex_normals(mesh::meshClass)
    I = vec(repeat(mesh.faces, 3, 1))
    println("Type of I: ", typeof(I))
    println("Shape of I: ", size(I))
    println("Value of I: ", I[3,1])
    J = vec(repeat(1:3, inner=3 * mesh.nf))
    println("Type of J: ", typeof(J))
    println("Shape of J: ", size(J))
    println("Value of J: ", J[43621,1])
    TA = spdiagm(repeat(mesh.ta, 3))
    println("Type of TA: ", typeof(TA))
    println("Shape of TA: ", size(TA))
    println("Value of TA: ", TA[3,2])
    println("Type of ta: ", typeof(mesh.ta))
    println("Shape of ta: ", size(mesh.ta))
    println("Value of ta: ", mesh.ta[3,1])
    S = vec(repeat(TA[1:mesh.nf, 1:mesh.nf] * mesh.Nf, 3))
    println("Type of S: ", typeof(S))
    println("Shape of S: ", size(S))
    println("Value of S: ", S[3,1])
    # Nv = sparse(I, J, S, mesh.nv, 3)
    # Nv = Array(sparse(I, J, S, mesh.nv, 3))
    # Nv = sparse(I, J, S, mesh.nv, 3)
    Nv = SparseArrays.sparse!(I, J, S, mesh.nv, 3)
    Nv = Matrix(Nv)
    println("Type of NV: ", typeof(Nv))
    println("Shape of NV: ", size(Nv))
    println("Value of NV: ", Nv[2,2])
    return normalize_vf(Nv)
end

function nc_data(mesh::meshClass)
    T = convert(Matrix{Int64}, mesh.faces)
    I = vec(hcat(T[:, 2], T[:, 3], T[:, 1]))
    J = vec(hcat(T[:, 3], T[:, 1], T[:, 2]))
    S = repeat(1:mesh.nf, inner=3)
    E = sparse(I, J, S, mesh.nv, mesh.nv)

    Elisto = hcat(I, J)
    sElist = mapslices(sort!, Elisto, dims=2)
    s = normv(Elisto - sElist) .> 1e-12
    t = S' .* (-1) .^ s
    # edges, une, _ = unique(sElist, dims=1, return_indices=true)
    edges = unique(sElist, dims=1)
    # _, une = findall(x -> x in edges, eachrow(sElist))
    _, une = findall(x -> x ∈ collect(eachrow(edges)), eachrow(sElist))
    ne = size(edges, 1)
    e2t = zeros(Int64, ne, 4)
    t2e = zeros(Int64, mesh.nf, 3)
    ie = zeros(Int64, ne)
    for m in 1:length(edges)
        i = edges[m, 1]
        j = edges[m, 2]
        t1 = t[une[m]]
        t2 = -(E[i, j] + E[j, i] - abs(t1)) * sign(t1)
        e2t[m, 1:2] = [t1, t2]
        f = T[abs(t1), :]
        loc = findfirst(x -> x == (sum(f) - i - j), f)
        t2e[abs(t1), loc] = m * sign(t1)
        e2t[m, 3] = loc
        if t2 != 0
            f = T[abs(t2), :]
            loc = findfirst(x -> x == (sum(f) - i - j), f)
            t2e[abs(t2), loc] = m * sign(t2)
            e2t[m, 4] = loc
            ie[m] = 1
        end
    end

    v2e = sparse(edges[:, 1], edges[:, 2], 1:length(edges), mesh.nv, mesh.nv)

    ne = size(edges, 1)
    nie = sum(ie)
    inner_edges = findall(x -> x != 0, ie)
    bv = zeros(Int64, mesh.nv)
    bv[vec(edges[ie.==0, :])] .= 1
    bf = zeros(Int64, mesh.nf)
    bf[sum(ismember(mesh.faces, findall(x -> x == 1, bv)), dims=2).>0] .= 1

    t1 = abs.(e2t[inner_edges, 1])
    t2 = abs.(e2t[inner_edges, 2])

    I = 1:2*nie
    S = ones(Float64, 2 * nie)
    e2t1 = sparse(I, vcat(t1, t1 + mesh.nf), S, 2 * nie, 2 * mesh.nf)
    e2t2 = sparse(I, vcat(t2, t2 + mesh.nf), S, 2 * nie, 2 * mesh.nf)
    return edges, e2t, t2e, e2t1, e2t2, v2e, ie, ne, nie, inner_edges, bv, bf
end

function compute_LB!(mesh::meshClass)
    mesh.Ww = cotLaplacian(mesh)
    laplacian = spdiagm(1 ./ mesh.va) * mesh.Ww
    mesh.Lap = laplacian
    mesh.Aa = Diagonal(mesh.va)
end

function GG!(mesh::meshClass)
    I = vec(repeat(1:mesh.nf, inner=3))
    II = vcat(I, I + mesh.nf, I + 2 * mesh.nf)
    J = vec(mesh.faces')
    JJ = vcat(J, J, J)
    RE1 = rotate_vf(mesh, mesh.E1)
    RE2 = rotate_vf(mesh, mesh.E2)
    RE3 = rotate_vf(mesh, mesh.E3)
    TA = mesh.ta
    S = vec(vcat(-RE1, RE2, -RE3))

    G = sparse(II, JJ, S, 3 * mesh.nf, mesh.nv)
    ITA = spdiagm(0.5 * vec(repeat(1 ./ TA, inner=3)))

    grad_op = ITA * G
    if any(isnan, grad_op)
        warn("Grad: NANs exist")
    end
    grad_op[isnan.(grad_op)] .= 0
    mesh.G = grad_op
end

function D!(mesh::meshClass)
    IVA = spdiagm(1 ./ mesh.va)
    TA = spdiagm(vec([mesh.ta; mesh.ta; mesh.ta]))
    D = -IVA * mesh.G' * TA
    mesh.D = D
end

function rotate_vf(mesh::meshClass, vf)
    vf = reshape(vf, mesh.nf, 3)
    return cross(mesh.Nf, vf)
end

# Function to calculate normal vectors norm
function normv(vf)
    return sqrt.(sum(vf .^ 2, dims=2))
end

# Function to normalize vector fields
# function normalize_vf(vf)
#     norm = normv(vf)
#     norm[norm.<1e-15] .= 1
#     return vf ./ norm
# end
function normalize_vf(vf)
    norm_vf = normv(vf)
    nnv = vf ./ hcat(norm_vf, norm_vf, norm_vf)
    nnv[norm_vf .< 1e-15, :] .= 0
    return nnv
end

# Function to get the camera settings
function get_camera(ca=gca())
    cam = Dict(
        :pba => get(ca, "PlotBoxAspectRatio"),
        :dar => get(ca, "DataAspectRatio"),
        :cva => get(ca, "CameraViewAngle"),
        :cuv => get(ca, "CameraUpVector"),
        :ct => get(ca, "CameraTarget"),
        :cp => get(ca, "CameraPosition")
    )
    return cam
end

# Function to set the camera settings
function set_camera(ca, cam)
    set(ca, "PlotBoxAspectRatio", cam[:pba])
    set(ca, "DataAspectRatio", cam[:dar])
    set(ca, "CameraViewAngle", cam[:cva])
    set(ca, "CameraUpVector", cam[:cuv])
    set(ca, "CameraTarget", cam[:ct])
    set(ca, "CameraPosition", cam[:cp])
end

end # module MeshClass
