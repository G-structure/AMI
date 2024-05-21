import numpy as np
from scipy.sparse import diags, spdiags, csr_matrix
from scipy.sparse.linalg import spsolve
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
            'H': Dirichlet (same as 'D', but with a different weight calculation)
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
    elif reg == 'H':
        alpha = alpha_hat * np.sqrt(np.sum(va) ** 3)
        Ww_s = curved_hessian(vertices, faces)  # Ensure that the mex `curved_hessian` is available
        varRho = 0
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
    QUIET = 1
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

    u_p = np.zeros(nv - len(x0))
    y = np.zeros(3*nf)
    z = np.zeros(3*nf)
    div_y = np.zeros(nv - len(x0))
    div_z = np.zeros(nv - len(x0))

    history = {
        'r_norm': np.zeros(niter),
        's_norm': np.zeros(niter),
        'eps_pri': np.zeros(niter),
        'eps_dual': np.zeros(niter)
    }

    # Eliminating x0 (b.c)
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

    if reg == 'vfa' or reg == 'H':
        Ww_s_p = Ww_s.copy()
        Ww_s_p[:, x0] = 0
        Ww_s_p[x0, :] = 0

    if not QUIET:
        print('{:3s}\t{:10s}\t{:10s}\t{:10s}\t{:10s}'.format(
            'iter', 'r norm', 'eps pri', 's norm', 'eps dual'))

    # Pre-factorization
    if reg == 'D':
        L = np.linalg.cholesky(Ww_p.toarray())
        P = np.arange(nv - len(x0))
    else:  # 'H', 'vfa'
        if not varRho:
            L = np.linalg.cholesky((alpha*Ww_s_p + rho*Ww_p).toarray())
            P = np.arange(nv - len(x0))

    for ii in range(niter):
        # Step 1 - u-minimization
        b = va_p - div_y + rho*div_z

        if reg == 'D':
            u_p = np.linalg.solve(L.T, np.linalg.solve(L, b[P]))[P] / (alpha + rho)
        else:  # 'H', 'vfa'
            if not varRho:
                u_p = np.linalg.solve(L.T, np.linalg.solve(L, b[P]))[P]
            else:  # reg == 'H' and varRho
                u_p = spsolve(alpha*Ww_s_p + rho*Ww_p, b)
        Gx = G_p @ u_p

        # Step 2 - z-minimization
        zold = z.copy()
        div_zold = div_z.copy()
        z = (1/rho)*y + Gx
        z = z.reshape(nf, 3)
        norm_z = np.sqrt(np.sum(z**2, axis=1))
        norm_z[norm_z < 1] = 1
        z = (z.T / norm_z).T
        z = z.ravel()
        div_z = div_p @ z

        # Step 3 - dual variable update
        y = y + rho*(alphak*Gx + (1-alphak)*zold - z)
        div_y = div_p @ y

        # Residuals update
        tasqGx = tasq * Gx
        tasqZ = tasq * z
        history['r_norm'][ii] = np.linalg.norm(tasqGx - tasqZ)
        history['s_norm'][ii] = rho * np.linalg.norm(div_z - div_zold)
        history['eps_pri'][ii] = thresh1 + RELTOL * max(np.linalg.norm(tasqGx), np.linalg.norm(tasqZ))
        history['eps_dual'][ii] = thresh2 + RELTOL * np.linalg.norm(div_y)

        if not QUIET:
            print('{:3d}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}'.format(
                ii, history['r_norm'][ii], history['eps_pri'][ii],
                history['s_norm'][ii], history['eps_dual'][ii]))

        if history['r_norm'][ii] < history['eps_pri'][ii] and history['s_norm'][ii] < history['eps_dual'][ii]:
            break

        # Varying penalty parameter
        if varRho:
            if history['r_norm'][ii] > mu * history['s_norm'][ii]:
                rho *= tauinc
                y /= tauinc
            elif history['s_norm'][ii] > mu * history['r_norm'][ii]:
                rho /= taudec
                y *= taudec

    u = np.zeros(nv)
    u[nv_p] = u_p

    return u