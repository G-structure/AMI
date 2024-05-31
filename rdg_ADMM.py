import numpy as np
from scipy.sparse import diags, spdiags, csr_matrix
from scipy.linalg import cho_factor, cho_solve
from MeshClass import MeshClass
from smoothVF import smooth_vf

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
    tasq = np.repeat(np.sqrt(ta), 3)

    # Set parameters according to the desired regularizer
    if reg == 'D':
        alpha = alpha_hat * np.sqrt(np.sum(va))
        varRho = 1
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
    niter = 10000
    QUIET = 0
    ABSTOL = 1e-5 / 2
    RELTOL = 1e-2
    mu = 10
    tauinc = 2
    taudec = 2
    alphak = 1.7

    if reg == 'H':
        ABSTOL /= 20
        RELTOL /= 20

    thresh1 = np.sqrt(3*nf) * ABSTOL * np.sqrt(np.sum(va))
    thresh2 = np.sqrt(nv) * ABSTOL * np.sum(va)

    nv_p = np.arange(nv)
    nv_p = np.delete(nv_p, x0)
    va_p = va.copy()
    va_p[x0] = 0
    Ww_p = Ww.copy()
    Ww_p[:, x0] = 0
    Ww_p[x0, :] = 0
    G_p = G.copy()
    G_p[:, x0] = 0
    G_pt = G_p.T
    div_p = G_pt.multiply(np.repeat(ta, 3))

    u_p = np.zeros(nv - len(x0))
    y = np.zeros(3*nf)
    z = np.zeros(3*nf)
    div_y = np.zeros(nv)
    div_z = np.zeros(nv)

    if reg == 'vfa' or reg == 'H':
        Ww_s_p = Ww_s.copy()
        Ww_s_p[:, x0] = 0
        Ww_s_p[x0, :] = 0

    if not QUIET:
        print('{:3s}\t{:10s}\t{:10s}\t{:10s}\t{:10s}'.format(
            'iter', 'r norm', 'eps pri', 's norm', 'eps dual'))

    # Pre-factorization
    if reg == 'D':
        Ww_p = Ww_p + 1e-10 * diags(np.ones(Ww_p.shape[0]))  # Add small regularization
        if not np.allclose(Ww_p.toarray(), Ww_p.toarray().T):
            raise ValueError("Ww_p is not symmetric")
        c, low = cho_factor(Ww_p.toarray())
    else:  # 'H', 'vfa'
        if not varRho:
            c, low = cho_factor(alpha * Ww_s_p + rho * Ww_p)

    for ii in range(niter):
        # Step 1 - u-minimization
        b = va_p - div_y + rho * div_z

        if reg == 'D':
            u_p = cho_solve((c, low), b) / (alpha + rho)
        else:  # 'H', 'vfa'
            if not varRho:
                u_p = cho_solve((c, low), b)
            else:  # 'H' with varRho
                u_p = np.linalg.solve(alpha * Ww_s_p + rho * Ww_p, b)

        Gx = G_p @ u_p

        # Step 2 - z-minimization
        zold = z.copy()
        div_zold = div_z.copy()
        z = (1 / rho) * y + Gx
        z = z.reshape(nf, 3).T
        max_abs_z = np.max(np.abs(z), axis=0)
        norm_z = max_abs_z * np.linalg.norm(z / max_abs_z, axis=0)
        norm_z[norm_z < 1] = 1
        z = (z / norm_z).T.ravel()
        div_z = div_p @ z

        # Step 3 - dual variable update
        y = y + rho * (alphak * Gx + (1 - alphak) * zold - z)
        div_y = div_p @ y

        # Residuals update
        tasqGx = tasq * Gx
        tasqZ = tasq * z
        max_abs_diff = np.max(np.abs(tasqGx - tasqZ))
        history_r_norm = max_abs_diff * np.linalg.norm(((tasqGx - tasqZ) / max_abs_diff).reshape(-1, 1), 'fro')

        max_abs_div_diff = np.max(np.abs(div_z - div_zold))
        history_s_norm = rho * max_abs_div_diff * np.linalg.norm(((div_z - div_zold) / max_abs_div_diff).reshape(-1, 1), 'fro')

        max_abs_tasqGx = np.max(np.abs(tasqGx))
        max_abs_tasqZ = np.max(np.abs(tasqZ))
        history_eps_pri = thresh1 + RELTOL * max(
            max_abs_tasqGx * np.linalg.norm((tasqGx / max_abs_tasqGx).reshape(-1, 1), 'fro'),
            max_abs_tasqZ * np.linalg.norm((tasqZ / max_abs_tasqZ).reshape(-1, 1), 'fro')
        )

        max_abs_div_y = np.max(np.abs(div_y))
        history_eps_dual = thresh2 + RELTOL * max_abs_div_y * np.linalg.norm((div_y / max_abs_div_y).reshape(-1, 1), 'fro')


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