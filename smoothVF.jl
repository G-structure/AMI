module smoothVF

using SparseArrays

export smooth_vf, ff

function smooth_vf(Mm, vf, n)
    # n - power vector; 2 for line fields

    # Based on code written by
    # "PH-CPF: Planar Hexagonal Meshing using Coordinate Power Fields" by [Pluta et al., 2021]

    nf = Mm.nf
    we = []

    locs = findall(x -> x > 1e-5, MeshClass.normv(vf))
    nl = length(locs)

    Aeq = sparse(collect(1:2*nl), vcat(locs, nf .+ locs), ones(2 * nl), 2 * nl, 2 * nf)

    # -> local basis -> power n -> locs
    beq = reshape(ff(Mm.EB * vf[:], n), :, 2)
    beq = beq[locs, :]
    beq = beq[:]

    C = Mm.godf(n)
    vf = zeros(2 * nf)

    # If there are no constraints; use the eigenvector
    if size(Aeq, 1) > 0
        x = lsqlin(C, vf, [], [], Aeq, beq)
    else
        x, _ = eigs(C, 1, "SM")
    end

    # -> sqrt n -> extrinsic
    w = reshape(Mm.EBI * ff(x, 1 / n), :, 3)
    # Normalize
    w = MeshClass.normalize_vf(w)

    return w
end

# in: 2n x 1
# out: 2n x 1, f[x,n] = x.^n
function ff(in, n)
    s = length(in) ÷ 2
    a = in[1:s]
    b = in[s+1:end]
    c = a .+ 1im .* b
    cn = c .^ n
    e = [real(cn); imag(cn)]
    return e
end

end # module smoothVF
