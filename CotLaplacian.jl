using SparseArrays

module CotLaplacian

export cotLaplacian

# Compute the cotangent weight Laplacian.
# W is the symmetric cot Laplacian; & A are the area weights
function cotLaplacian(mesh, L23, L13, L12)
    X = mesh.vertices
    T = mesh.faces
    nv = size(X,1)

    inputL = (varinfo(cotLaplacian).name .== "L23") +
             (varinfo(cotLaplacian).name .== "L13") +
             (varinfo(cotLaplacian).name .== "L12")

    if inputL .< 3
        # Find orig edge lengths & angles
        L1 = normv(X[T[:,2],:] .- X[T[:,3],:])
        L2 = normv(X[T[:,1],:] .- X[T[:,3],:])
        L3 = normv(X[T[:,1],:] .- X[T[:,2],:])
    else
        L1 = L23
        L2 = L13
        L3 = L12
    end

    EL = [L1, L2, L3]
    A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2 .* L2 .* L3)
    A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2 .* L1 .* L3)
    A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2 .* L1 .* L2)
    A = [A1, A2, A3]
    A = acos.(A)

    # The Cot Laplacian
    I = [T[:,1] T[:,2] T[:,3]]'
    J = [T[:,2] T[:,3] T[:,1]]'
    S = 0.5 * cot.(collect([A[:,3] A[:,1] A[:,2]]'))
    In = vcat(I, J, I, J)'
    Jn = vcat(J, I, I, J)'
    Sn = vcat(-S, -S, S, S)'
    W = sparse(In, Jn, Sn, nv, nv)

    if inputL .< 3
        # Use the Barycentric areas
        M = mass_matrix_barycentric(mesh)
        A = sum(M, dims=2)
    else
        M = mass_matrix_barycentric(mesh, L1, L2, L3)
        A = sum(M, dims=2)
    end

    return W, A
end

function normv(V)
    sqrt.(sum(V.^2, dims=2))
end

function mass_matrix_barycentric(mesh, L1=nothing, L2=nothing, L3=nothing)
    T = mesh.faces;
    inputL = (L1 == nothing) + (L2 == nothing) + (L3 == nothing)

    if inputL >= 3
        Ar = mesh.ta
    else
        s = (L1 + L2 + L3) / 2
        Ar = sqrt.(s .* (s .- L1) .* (s .- L2) .* (s .- L3));
    end

    nv = mesh.nv

    I = [T[:,1] T[:,2] T[:,3]]'
    J = [T[:,2] T[:,3] T[:,1]]'
    Mij = (1 / 12) * vcat(Ar, Ar, Ar)'
    Mji = Mij
    Mii = (1 / 6) * vcat(Ar, Ar, Ar)'
    In = vcat(I, J, I)'
    Jn = vcat(J, I, I)'
    Mn = vcat(Mij, Mji, Mii)'
    M = sparse(In, Jn, Mn, nv, nv)

    return M
end

end  # module cotLaplacian
