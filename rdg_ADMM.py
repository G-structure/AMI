import numpy as np
from scipy.sparse import diags, spdiags, csr_matrix, csc_matrix
from scipy.linalg import cho_factor, cho_solve
from MeshClass import MeshClass
from smoothVF import smooth_vf
from scipy.sparse.csgraph import reverse_cuthill_mckee

def rdg_ADMM(Mm, x0, reg='D', alpha_hat=0.1, beta_hat=0, vf=0):
    """
    ADMM algorithm for computing regularized geodesic distances.

    Parameters:
        Mm (MeshClass): An instance of MeshClass containing the mesh data.
        x0 (list or numpy.ndarray): A vector of source set vertex indices.
        reg (str, optional): The type of regularizer to use. Options are:
            'D': Dirichlet (default)
            'vfa': Vector Field Alignment
        alpha_hat (float, optional): The regularizer weight, which is scale-invariant. Default is 0.1.
        beta_hat (float, optional): The vector field alignment weight, relevant when `reg` is 'vfa'. Default is 0.
        vf (numpy.ndarray, optional): A matrix of size |F|x3 representing the vector field to align to, relevant when `reg` is 'vfa'. Default is 0.

    Returns:
        u (numpy.ndarray): The computed regularized distance.
    """
    # Mesh data
    vertices = Mm.vertices
    faces = Mm.faces
    nv = Mm.nv
    nf = Mm.nf
    va = Mm.va
    ta = Mm.ta
    G = Mm.G
    Ww = Mm.Ww

      # Mesh data
    vertices = Mm.vertices
    faces = Mm.faces
    nv = Mm.nv
    nf = Mm.nf
    va = Mm.va
    ta = Mm.ta
    G = Mm.G
    Ww = Mm.Ww

    print(f"vertices shape: {vertices.shape}")
    print(f"faces shape: {faces.shape}")
    print(faces)
    print(f"nv shape: {nv if isinstance(nv, int) else nv.shape}")
    print(f"nf shape: {nf if isinstance(nf, int) else nf.shape}")
    print(f"va shape: {va.shape}")
    print(f"ta shape: {ta.shape}")
    print(f"G shape: {G.shape}")
    print(f"Ww shape: {Ww.shape}")

    # Log the 5th element for each column of vertices and faces
    if vertices.shape[0] > 4:
        for i in range(vertices.shape[1]):
            print(f"vertices[4, {i}] value: {vertices[4, i]}")
    if faces.shape[0] > 4:
        for i in range(faces.shape[1]):
            print(f"faces[4, {i}] value: {faces[4, i]}")

    # Log the 5th element for va, ta, tasq
    if va.shape[0] > 4:
        print(f"va[4] value: {va[4]}")
    if ta.shape[0] > 4:
        print(f"ta[4] value: {ta[4]}")
    tasq = np.repeat(np.sqrt(ta), 3)
    print(f"tasq shape: {tasq.shape}")
    if tasq.shape[0] > 4:
        print(f"tasq[4] value: {tasq[4]}")

    # Log the 5th element of the first column for G and Ww
    if G.shape[0] > 4:
        print(f"G[4, 0] value: {G[4, 0]}")
    if Ww.shape[0] > 4:
        print(f"Ww[4, 0] value: {Ww[4, 0]}")

    # Set parameters according to the desired regularizer
    if reg == 'D':
        alpha = alpha_hat * np.sqrt(np.sum(va))
        varRho = 1
        print(f"alpha shape: {alpha.shape if hasattr(alpha, 'shape') else 'scalar'}")
        print(f"alpha value: {alpha}")
    elif reg == 'vfa':
        alpha = alpha_hat * np.sqrt(np.sum(va))
        beta = beta_hat * np.sqrt(np.sum(va))
        if np.max(MeshClass.normv(vf)) < 1e-10:
            raise ValueError("Vector field for alignment is empty")
        Vmat = np.vstack((
            spdiags(vf[:, 0] * vf[:, 0], 0, nf, nf),
            spdiags(vf[:, 0] * vf[:, 1], 0, nf, nf),
            spdiags(vf[:, 0] * vf[:, 2], 0, nf, nf),
            spdiags(vf[:, 1] * vf[:, 0], 0, nf, nf),
            spdiags(vf[:, 1] * vf[:, 1], 0, nf, nf),
            spdiags(vf[:, 1] * vf[:, 2], 0, nf, nf),
            spdiags(vf[:, 2] * vf[:, 0], 0, nf, nf),
            spdiags(vf[:, 2] * vf[:, 1], 0, nf, nf),
            spdiags(vf[:, 2] * vf[:, 2], 0, nf, nf)
        ))
        Ww_s = G.T @ spdiags(np.repeat(ta, 3), 0, 3*nf, 3*nf) @ (csr_matrix(np.eye(3*nf)) + beta*Vmat) @ G
        varRho = 0
    else:
        raise ValueError("Unrecognized regularizer")

    # ADMM parameters
    rho = 2 * np.sqrt(np.sum(va))
    print(rho)
    niter = 100
    QUIET = 0
    ABSTOL = 1e-5 / 2
    RELTOL = 1e-2
    mu = 10
    tauinc = 2
    taudec = 2
    alphak = 1.7

    thresh1 = np.sqrt(3*nf) * ABSTOL * np.sqrt(np.sum(va))
    thresh2 = np.sqrt(nv) * ABSTOL * np.sum(va)

    # nv_p = np.arange(nv)
    # nv_p = np.delete(nv_p, x0)
    # va_p = va.copy()
    # va_p[x0] = 0
    # Ww_p = Ww.copy()
    # Ww_p[:, x0] = 0
    # Ww_p[x0, :] = 0
    # G_p = G.copy()
    # G_p[:, x0] = 0
    # G_pt = G_p.T
    # div_p = G_pt.multiply(np.repeat(ta, 3))

    u_p = np.zeros((nv - len(x0), 1))
    y = np.zeros((3*nf, 1))
    z = np.zeros((3*nf, 1))
    div_y = np.zeros((nv - len(x0), 1))
    div_z = np.zeros((nv - len(x0), 1))

    print("u_p shape:", u_p.shape)

    nv_p = np.arange(nv)
    nv_p = np.delete(nv_p, x0)
    va_p = va[nv_p]
    Ww_p = Ww[nv_p, :][:, nv_p]
    G_p = G[:, nv_p]
    G_pt = G_p.T
    div_p = G_pt.multiply(np.repeat(ta, 3))

    if reg == 'vfa':
        Ww_s_p = Ww_s.copy()
        Ww_s_p[:, x0] = 0
        Ww_s_p[x0, :] = 0

    if not QUIET:
        print('{:3s}\t{:10s}\t{:10s}\t{:10s}\t{:10s}'.format(
            'iter', 'r norm', 'eps pri', 's norm', 'eps dual'))

## OLD CODE
    # # Pre-factorization
    # if reg == 'D':
    #     Ww_p = Ww_p + 1e-10 * diags(np.ones(Ww_p.shape[0]))  # Add small regularization
    #     print(f"Ww_p shape: {Ww_p.shape}")

    #     if not np.allclose(Ww_p.toarray(), Ww_p.toarray().T):
    #         raise ValueError("Ww_p is not symmetric")

    #     c, low = cho_factor(Ww_p.toarray())
    #     print(f"c shape: {c.shape}")
    #     print(f"low value: {low}")
    # else:  # 'H', 'vfa'
    #     if not varRho:
    #         print(f"alpha shape: {alpha.shape}")
    #         print(f"Ww_s_p shape: {Ww_s_p.shape}")
    #         print(f"rho shape: {rho.shape}")
    #         print(f"Ww_p shape: {Ww_p.shape}")

    #         c, low = cho_factor(alpha * Ww_s_p + rho * Ww_p)
    #         print(f"c shape: {c.shape}")
    #         print(f"low value: {low}")

    # for ii in range(niter):
    #     # Step 1 - u-minimization
    #     b = va_p - div_y + rho * div_z

    #     if reg == 'D':
    #         u_p = cho_solve((c, low), b) / (alpha + rho)
    #     else:  # 'H', 'vfa'
    #         if not varRho:
    #             u_p = cho_solve((c, low), b)
    #         else:  # 'H' with varRho
    #             u_p = np.linalg.solve(alpha * Ww_s_p + rho * Ww_p, b)

    #     Gx = G_p @ u_p




    # Pre-factorization
    if reg == 'D':
        Ww_p_dense = Ww_p.toarray()
        perm = reverse_cuthill_mckee(csc_matrix(Ww_p_dense))
        P = np.eye(Ww_p_dense.shape[0])[:, perm]
        Ww_p_permuted = P.T @ Ww_p_dense @ P
        c, low = cho_factor(Ww_p_permuted, lower=True)
        print(f"Size of L: {c.shape}")
        print(f"Size of P: {P.shape}")
    else:  # 'H', 'vfa'
        if not varRho:
            Ww_combined = alpha * Ww_s_p + rho * Ww_p
            Ww_combined_dense = Ww_combined.toarray()
            perm = reverse_cuthill_mckee(csc_matrix(Ww_combined_dense))
            P = np.eye(Ww_combined_dense.shape[0])[:, perm]
            Ww_combined_permuted = P.T @ Ww_combined_dense @ P
            c, low = cho_factor(Ww_combined_permuted, lower=True)
            print(f"Size of L: {c.shape}")
            print(f"Size of P: {P.shape}")

    for ii in range(niter):
        # Step 1 - u-minimization
        print("div_y shape:", div_y.shape)
        print("div_z shape:", div_z.shape)
        print("div_p shape:", div_p.shape)
        print("va_p shape:", va_p.shape)
        b = va_p[:, np.newaxis] - div_y + rho * div_z
        print("b shape:", b.shape)


        if reg == 'D':
            u_p = P @ cho_solve((c, low), P.T @ b) / (alpha + rho)
        else:  # 'H', 'vfa'
            if not varRho:
                u_p = P @ cho_solve((c, low), P.T @ b)
            else:  # 'H' with varRho
                u_p = np.linalg.solve(alpha * Ww_s_p + rho * Ww_p, b)


        Gx = G_p @ u_p


        # Step 2 - z-minimization
        # print("y shape:", y.shape)
        # print("Gx shape:", Gx.shape)
        # print("Gp shape:", G_p.shape)
        # print("u_p shape:", u_p.shape)
        # print("div_p shape:", div_p.shape)
        # print("div_y shape:", div_y.shape)
        # print("div_z shape:", div_z.shape)
        # print("va_p shape:", va_p.shape)
        # print("z shape:", z.shape)
        # print("rho shape:", rho.shape)
        zold = z.copy()
        div_zold = div_z.copy()
        z = (1 / rho) * y + Gx
        print("z shape:", z.shape)
        # print(nf)
        z = np.reshape(z, (nf, 3), order='F').T
        # z = z.reshape((3, nf), order='F').T
        print("z shape:", z.shape)
        norm_z = np.sqrt(np.sum(z**2, axis=0, keepdims=True))
        print("norm_z shape:", norm_z.shape)
        norm_z[norm_z < 1] = 1
        print("norm_z shape:", norm_z.shape)
        # z = (z / norm_z).T.ravel()
        z = (z / norm_z).T.ravel(order='F')
        print("z shape:", z.shape)
        div_z = (div_p @ z).reshape(-1, 1)
        print("div_z shape:", div_z.shape)

        # Step 3 - dual variable update
        y = y + rho * (alphak * Gx + (1 - alphak) * zold - z)
        div_y = div_p @ y

        # Residuals update
        tasqGx = tasq * Gx
        tasqZ = tasq * z
        print("hi starting long funct")
        history_r_norm = np.linalg.norm(tasqGx - tasqZ)
        print("r norm done")
        history_s_norm = rho * np.linalg.norm(div_z - div_zold)
        print("s norm done")
        history_eps_pri = thresh1 + RELTOL * max(np.linalg.norm(tasqGx), np.linalg.norm(tasqZ))
        print("eps pri done")
        history_eps_dual = thresh2 + RELTOL * np.linalg.norm(div_y)
        print("eps dual done")

        if not QUIET:
            print(f'{ii:3d}\t{history_r_norm:10.4f}\t{history_eps_pri:10.4f}\t{history_s_norm:10.4f}\t{history_eps_dual:10.4f}')

        # Stopping criteria
        if ii > 1 and history_r_norm < history_eps_pri and history_s_norm < history_eps_dual:
            break

        # Varying penalty parameter
        if varRho:
            if history_r_norm > mu * history_s_norm:
                rho *= tauinc
            elif history_s_norm > mu * history_r_norm:
                rho /= taudec

    u = np.zeros(nv)
    u[nv_p] = u_p

    return u
