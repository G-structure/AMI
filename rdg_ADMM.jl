module RDG_ADMM

using SparseArrays
using LinearAlgebra

include("MeshClass.jl")
using .MeshClass

export rdg_ADMM

function rdg_ADMM(Mm::meshClass, x0; reg="D", alpha_hat=0.1, beta_hat=0.0, vf=nothing)
    # ADMM algorithm for computing regularized geodesic distances
    #
    # Inputs:
    #   Mm - MeshClass
    #   x0 - source set - vertex indices
    #   Optional:
    #       reg - which regularizer to use.
    #           'D' - dirichlet
    #           "vfa" - vector field alignment
    #       alpha_hat - regularizer weight; scale invariant
    #       beta_hat - vector field alignment weight; scale-invariant; relevant for reg = "vfa"
    #       vf - |F|x3 - vector field to align to; relevant for reg = "vfa"
    #
    # Outputs:
    #   u - the computed regularized distance

    reg = lowercase(reg)
    nv = Mm.nv
    nf = Mm.nf
    va = Mm.va
    ta = Mm.ta
    G = Mm.G
    Ww = Mm.Ww
    tasq = repeat(sqrt.(ta), 3, 1)

    # Set parameters according to the desired regularizer
    if reg == 'd'
        alpha = alpha_hat * sqrt(sum(va))
        varRho = true  # determine whether to use a varying penalty parameter
    elseif reg == "vfa"
        alpha = alpha_hat * sqrt(sum(va))
        beta = beta_hat * sqrt(sum(va))
        if maximum(norm.(vf, dims=2)) < 1e-10
            error("Vector field for alignment is empty")
        end
        Vmat = vcat(
            hcat(spdiagm(vf[:,1] .* vf[:,1]), spdiagm(vf[:,1] .* vf[:,2]), spdiagm(vf[:,1] .* vf[:,3])),
            hcat(spdiagm(vf[:,2] .* vf[:,1]), spdiagm(vf[:,2] .* vf[:,2]), spdiagm(vf[:,2] .* vf[:,3])),
            hcat(spdiagm(vf[:,3] .* vf[:,1]), spdiagm(vf[:,3] .* vf[:,2]), spdiagm(vf[:,3] .* vf[:,3]))
        )
        Ww_s = G' * spdiagm(repeat(ta, 3, 1)) * (I + beta * Vmat) * G
        varRho = false  # determine whether to use a varying penalty parameter
    else
        error("Unrecognized regularizer")
    end

    # ADMM parameters
    rho = 2 * sqrt(sum(va))
    niter = 10000
    QUIET = true
    ABSTOL = 1e-5 / 2
    RELTOL = 1e-2
    mu = 10.0  # >1
    tauinc = 2.0  # >1
    taudec = 2.0  # >1
    alphak = 1.7  # over-relaxation

    thresh1 = sqrt(3 * nf) * ABSTOL * sqrt(sum(va))
    thresh2 = sqrt(nv) * ABSTOL * sum(va)

    u_p = zeros(Float64, nv - length(x0))
    y = zeros(Float64, 3 * nf)
    z = zeros(Float64, 3 * nf)
    div_y = zeros(Float64, nv - length(x0))
    div_z = zeros(Float64, nv - length(x0))

    history_r_norm = zeros(Float64, niter)
    history_s_norm = zeros(Float64, niter)
    history_eps_pri = zeros(Float64, niter)
    history_eps_dual = zeros(Float64, niter)

    # Eliminating x0 (boundary conditions)
    nv_p = collect(1:nv)
    nv_p[x0] .= []
    va_p = va; va_p[x0] .= 0
    Ww_p = Ww; Ww_p[:, x0] .= 0; Ww_p[x0, :] .= 0
    G_p = G; G_p[:, x0] .= 0; G_pt = G_p'
    div_p = G_pt .* repeat(ta, 3, 1)'

    if reg == "vfa"
        Ww_s_p = Ww_s
        Ww_s_p[:, x0] .= 0
        Ww_s_p[x0, :] .= 0
    end

    if !QUIET
        println("| Iter |   r norm   |  eps pri  |   s norm   |  eps dual  |")
    end

    # Pre-factorization
    if reg == 'd'
        L, _, P = cholesky(Ww_p, :L)
    else
        if !varRho
            L, _, P = cholesky(alpha * Ww_s_p + rho * Ww_p, :L)
        end
    end

    for ii in 1:niter
        # Step 1 - u-minimization
        b = va_p .- div_y .+ rho .* div_z
        if reg == 'd'
            u_p = P * (L \ (L' \ (P' * b))) / (alpha + rho)
        else
            if !varRho
                u_p = P * (L \ (L' \ (P' * b)))
            end
        end
        Gx = G_p * u_p

        # Step 2 - z-minimization
        zold = z
        div_zold = div_z
        z = (1/rho) * y .+ Gx
        z = reshape(z, nf, 3)'
        norm_z = sqrt.(sum(z.^2, dims=2))
        norm_z[norm_z .< 1] .= 1
        z = z ./ norm_z
        z = z'; z = vec(z)
        div_z = div_p * z

        # Step 3 - dual variable update
        y = y .+ rho .* (alphak * Gx .+ (1 - alphak) * zold .- z)
        div_y = div_p * y

        # Residuals update
        tasqGx = tasq .* Gx
        tasqZ = tasq .* z
        history_r_norm[ii] = norm(tasqGx - tasqZ, "fro")
        history_s_norm[ii] = rho * norm(div_z - div_zold, "fro")
        history_eps_pri[ii] = thresh1 + RELTOL * max(norm(tasqGx, "fro"), norm(tasqZ, "fro"))
        history_eps_dual[ii] = thresh2 + RELTOL * norm(div_y, "fro")

        if !QUIET
            println("|", ii, "|", history_r_norm[ii], "|", history_eps_pri[ii], "|", history_s_norm[ii], "|", history_eps_dual[ii], "|")
        end

        # Stopping criteria
        if ii > 1 && (history_r_norm[ii] < history_eps_pri[ii] && history_s_norm[ii] < history_eps_dual[ii])
            break
        end

        # Varying penalty parameter
        if varRho
            if history_r_norm[ii] > mu * history_s_norm[ii]
                rho = tauinc * rho
            elseif history_s_norm[ii] > mu * history_r_norm[ii]
                rho = rho / taudec
            end
        end
    end

    u = zeros(Float64, nv)
    u[nv_p] = u_p

    return u
end

end # module rdg_ADMM
