import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import eigs
import MeshClass

def smooth_vf(Mm, vf, n):
    """
    Smooth a vector field on a mesh.

    Parameters:
        Mm (MeshClass): A mesh object containing mesh information.
        vf (numpy.ndarray): A vector field defined on the mesh that needs to be smoothed.
        n (float or numpy.ndarray): A scalar or vector that defines the power to which the local basis is raised during smoothing.

    Returns:
        w (numpy.ndarray): The smoothed vector field, normalized and reshaped to match the extrinsic dimensions of the mesh.
    """
    nf = Mm.nf
    we = []

    locs = np.where(MeshClass.normv(vf) > 1e-5)[0]
    nl = locs.shape[0]

    Aeq = csr_matrix((np.ones(2 * nl), (np.arange(2 * nl), np.concatenate((locs, nf + locs)))), shape=(2 * nl, 2 * nf))

    beq = ff(Mm.EB[locs] @ vf, n).reshape(-1, 2)
    beq = beq.ravel()

    C, _ = Mm.godf(n)
    vf = np.zeros(nf * 2)

    if Aeq.shape[0] > 0:
        x = lsq_linear(C, vf, Aeq, beq).x
    else:
        x = eigs(C, k=1, which='SM')[1]

    w = MeshClass.normalize_vf(Mm.EBI @ ff(x, 1 / n).reshape(-1, 3))

    return w


def ff(inp, n):
    """
    Element-wise operations on vectors.

    Parameters:
        inp (numpy.ndarray): A vector of complex numbers represented as a 2n x 1 real vector,
                             where the first half represents the real parts and the second half represents the imaginary parts.
        n (float): The power to which each complex number is raised.

    Returns:
        e (numpy.ndarray): A 2n x 1 real vector representing the real and imaginary parts of the complex numbers after being raised to the power `n`.
    """
    s = inp.shape[0] // 2
    a = inp[:s]
    b = inp[s:]
    c = a + 1j * b
    cn = c ** n
    e = np.concatenate((cn.real, cn.imag))
    return e