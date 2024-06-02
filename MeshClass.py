import numpy as np
from scipy.sparse import csr_matrix, diags, eye
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MeshClass:
    def __init__(self, vertices=None, faces=None, name=None):
        if vertices is not None and faces is not None:
            self.vertices = np.array(vertices, dtype=np.float64)
            self.faces = np.array(faces, dtype=np.int32)
            self.name = name if name is not None else "Unnamed Mesh"
        else:
            raise ValueError("Please provide vertices and faces.")
        self.compute_all()

    def compute_all(self):
        self.nv = self.vertices.shape[0]
        self.nf = self.faces.shape[0]
        Nf = np.cross(self.vertices[self.faces[:, 0]] - self.vertices[self.faces[:, 1]],
            self.vertices[self.faces[:, 0]] - self.vertices[self.faces[:, 2]])
        self.Nf = self.normalize_vf(Nf);
        self.ta = np.sqrt(np.sum(Nf**2, axis=1)) / 2
        self.Nv = self.vertex_normals()
        self.va = self.calculatefvConnectivity().T @ self.ta / 3
        self.E1 = self.vertices[self.faces[:, 1]] - self.vertices[self.faces[:, 2]]
        self.E2 = self.vertices[self.faces[:, 0]] - self.vertices[self.faces[:, 2]]
        self.E3 = self.vertices[self.faces[:, 0]] - self.vertices[self.faces[:, 1]]
        self.R = self.rot()
        EE = np.sort(np.concatenate((self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]])), axis=1)
        self.E, une = np.unique(EE, axis=0, return_index=True)
        self.edges, self.e2t, self.t2e, self.e2t1, self.e2t2, self.v2e, self.ie, self.ne, self.nie, self.inner_edges, self.bv, self.bf = self.nc_data()
        self.ea = self.edge_areas()
        self.edge_basis()
        self.compute_LB()
        self.GG()
        self.DD()

    def interpulateFace2Vertices(self, fF=None):
        Af = self.ta
        Av = self.va
        I_F2V = csr_matrix((np.concatenate((Af / 3 / Av[self.faces[:, 0]], Af / 3 / Av[self.faces[:, 1]], Af / 3 / Av[self.faces[:, 2]])),
                            (np.concatenate((np.arange(self.nf), np.arange(self.nf), np.arange(self.nf))),
                             np.concatenate((self.faces[:, 0], self.faces[:, 1], self.faces[:, 2])))),
                           shape=(self.nf, self.nv))
        if fF is not None:
            fv = I_F2V @ fF
            return fv, I_F2V
        else:
            return None, I_F2V

    def interpulateVertices2Face(self, fV=None):
        Afinv = diags(1 / self.ta)
        Av = diags(self.va)
        _, I_F2V = self.interpulateFace2Vertices(np.ones(self.nf))
        I_V2F = Afinv @ I_F2V.T @ Av
        if fV is not None:
            fF = I_V2F @ fV
            return fF, I_V2F
        else:
            return None, I_V2F

    def edge_basis(self):
        NE1 = self.normalize_vf(self.E1)
        NE2 = self.R @ NE1.reshape(-1, 1)
        NE2 = NE2.reshape(-1, 3)
        I = np.repeat(np.arange(self.nf), 3)
        J = np.tile(np.arange(self.nf), 3) + np.repeat(np.arange(3) * self.nf, self.nf)
        B1 = csr_matrix((NE1.ravel(), (I, J)), shape=(self.nf, 3 * self.nf))
        B2 = csr_matrix((NE2.ravel(), (I, J)), shape=(self.nf, 3 * self.nf))
        self.EB = csr_matrix(np.vstack((B1.toarray(), B2.toarray())))
        self.EBI = self.EB.T
        self.F1 = NE1
        self.F2 = NE2

    def rot(self):
        sf = self.nf
        n = self.Nf
        II = np.repeat(np.arange(sf), 6)
        JJ1 = np.arange(sf)
        JJ2 = JJ1 + sf
        JJ3 = JJ2 + sf
        JJ = np.concatenate((JJ2, JJ3, JJ1, JJ3, JJ1, JJ2))
        SS = np.concatenate((-n[:, 2], n[:, 1], n[:, 2], -n[:, 0], -n[:, 1], n[:, 0]))
        R = csr_matrix((SS, (II, JJ)), shape=(3 * sf, 3 * sf))
        return R

    def calculatefvConnectivity(self):
        fvConnectivity = csr_matrix((np.ones(3 * self.nf), (np.repeat(np.arange(self.nf), 3), self.faces.ravel())), shape=(self.nf, self.nv))
        return fvConnectivity

    def baryCentersCalc(self):
        v1 = self.vertices[self.faces[:, 0]]
        v2 = self.vertices[self.faces[:, 1]]
        v3 = self.vertices[self.faces[:, 2]]
        baryCenters = (v1 + v2 + v3) / 3
        return baryCenters

    def edge_areas(self):
        T = self.faces
        I = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))
        J = np.concatenate((T[:, 2], T[:, 0], T[:, 1]))
        S = np.repeat(self.ta / 3, 3)
        In = np.concatenate((I, J))
        Jn = np.concatenate((J, I))
        Sn = np.concatenate((S, S))
        W = csr_matrix((Sn, (In, Jn)), shape=(self.nv, self.nv))
        ea = np.zeros(self.edges.shape[0])
        for i in range(len(ea)):
            ea[i] = W[self.edges[i, 0], self.edges[i, 1]]
            ea[i] += W[self.edges[i, 1], self.edges[i, 0]]
        return ea

    def rotate_vf(self, vf):
        vf = vf.reshape(self.nf, 3)
        rvf = np.cross(self.Nf, vf)
        return rvf

    def compute_LB(self):
        self.Ww, _ = self.cotLaplacian(self)
        laplacian = diags(1 / self.va) @ self.Ww
        self.Lap = laplacian
        self.Aa = diags(self.va)

    def GG(self):
        I = np.tile(np.arange(self.nf), (3, 1))
        II = np.concatenate((I.ravel(), I.ravel() + self.nf, I.ravel() + 2 * self.nf))
        J = self.faces.T.ravel()
        JJ = np.repeat(J, 3)
        RE1 = self.rotate_vf(self.E1)
        RE2 = self.rotate_vf(self.E2)
        RE3 = self.rotate_vf(self.E3)
        TA = self.ta
        S = np.concatenate((-RE1.ravel(), RE2.ravel(), -RE3.ravel()))
        G = csr_matrix((S, (II, JJ)), shape=(3 * self.nf, self.nv))
        ITA = diags(0.5 / np.repeat(TA, 3))
        grad_op = ITA @ G
        if np.any(np.isnan(grad_op.data)):
            print("Grad: NANs exist")
        grad_op.data[np.isnan(grad_op.data)] = 0
        self.G = grad_op

    def DD(self):
        IVA = diags(1 / self.va)
        TA = diags(np.repeat(self.ta, 3))
        D = -IVA @ self.G.T @ TA
        self.D = D

    def vertex_normals(self):
        I = np.repeat(self.faces, 3, axis=0).ravel()
        J = np.tile(np.arange(3), 3 * self.nf)

        TA = diags(np.repeat(self.ta, 3))
        S = np.repeat(TA.diagonal()[:self.nf, np.newaxis] * self.Nf, 3, axis=0).ravel()
        Nv = csr_matrix((S, (I, J)), shape=(self.nv, 3)).toarray()
        Nv = self.normalize_vf(Nv)
        return Nv

    def nc_data(self):
        T = self.faces
        I = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))
        J = np.concatenate((T[:, 2], T[:, 0], T[:, 1]))
        S = np.arange(1, self.nf + 1)
        S = np.concatenate((S, S, S))
        E = csr_matrix((S, (I, J)), shape=(self.nv, self.nv))
        Elisto = np.vstack((I, J)).T
        sElist = np.sort(Elisto, axis=1)
        s = (self.normv(Elisto - sElist) > 1e-12)
        t = S * (-1) ** s
        edges, une = np.unique(sElist, axis=0, return_index=True)
        ne = edges.shape[0]
        e2t = np.zeros((ne, 4), dtype=int)
        t2e = np.zeros((self.nf, 3), dtype=int)
        ie = np.zeros(ne, dtype=int)
        for m in range(len(edges)):
            i, j = edges[m]
            t1 = t[une[m]]
            t2 = -(E[i, j] + E[j, i] - abs(t1)) * np.sign(t1)
            e2t[m, :2] = [t1, t2]
            f = T[abs(t1) - 1]
            loc = np.where(f == (f.sum() - i - j))[0]
            t2e[abs(t1) - 1, loc] = m * np.sign(t1)
            e2t[m, 2] = loc
            if t2 != 0:
                f = T[abs(t2) - 1]
                loc = np.where(f == (f.sum() - i - j))[0]
                t2e[abs(t2) - 1, loc] = m * np.sign(t2)
                e2t[m, 3] = loc
                ie[m] = 1
        v2e = csr_matrix((np.arange(1, len(edges) + 1), (edges[:, 0], edges[:, 1])), shape=(self.nv, self.nv))
        ne = edges.shape[0]
        nie = np.sum(ie)
        inner_edges = np.where(ie)[0]
        bv = np.zeros(self.nv, dtype=int)
        bv[edges[ie == 0].ravel()] = 1
        bf = np.zeros(self.nf, dtype=int)
        bf[np.sum(np.isin(self.faces, np.where(bv == 1)[0]), axis=1) > 0] = 1
        t1 = abs(e2t[inner_edges, 0])
        t2 = abs(e2t[inner_edges, 1])
        I = np.arange(2 * nie)
        S = np.ones(2 * nie)
        e2t1 = csr_matrix((S, (I, np.concatenate((t1 - 1, t1 + self.nf - 1)))), shape=(2 * nie, 2 * self.nf))
        e2t2 = csr_matrix((S, (I, np.concatenate((t2 - 1, t2 + self.nf - 1)))), shape=(2 * nie, 2 * self.nf))
        return edges, e2t, t2e, e2t1, e2t2, v2e, ie, ne, nie, inner_edges, bv, bf

    def normalize_mesh(self, bbdO=1):
        xx = self.vertices - self.vertices.mean(axis=0)
        bbd = np.linalg.norm(self.vertices.max(axis=0) - self.vertices.min(axis=0))
        xx = xx / bbd * bbdO
        self.vertices = xx
        self.compute_all()

    def center_mesh(self, aa=None):
        if aa is not None:
            xx = self.vertices - aa
        else:
            xx = self.vertices - self.vertices.mean(axis=0)
            aa = self.vertices.mean(axis=0)
        self.vertices = xx
        self.compute_all()
        return aa

    def scale_mesh(self, scale_fac):
        xx = self.vertices * scale_fac
        self.vertices = xx
        self.compute_all()

    def visualizeMesh(self, f_vertices=None, f_faces=None, edgeColorFlag=0, figFlag=1):
        if f_vertices is None:
            f_vertices = np.array([])
        if f_faces is None:
            f_faces = np.array([])
        f_vertices = f_vertices.ravel()
        f_faces = f_faces.ravel()
        if figFlag:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = plt.gca()
        if f_vertices.size > 0 and f_faces.size > 0:
            print("Cannot display f_vertices and f_faces")
        elif f_vertices.size > 0:
            p = ax.plot_trisurf(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2],
                                triangles=self.faces, cmap='viridis', edgecolor='none' if edgeColorFlag else 'k')
            p.set_array(f_vertices)
            ax.set_title(self.name)
            fig.colorbar(p)
        elif f_faces.size > 0:
            p = ax.plot_trisurf(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2],
                                triangles=self.faces, cmap='viridis', edgecolor='none' if edgeColorFlag else 'k')
            color = np.ones(self.nv) * np.nan
            color[self.faces.ravel()] = f_faces
            p.set_array(color)
            ax.set_title(self.name)
            fig.colorbar(p)
        else:
            p = ax.plot_trisurf(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2],
                                triangles=self.faces, color='w', alpha=0.5, edgecolor='none' if edgeColorFlag else 'k')
            ax.set_title(self.name)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.view_init(elev=90, azim=0)
        plt.tight_layout()
        return p

    def vectorFieldVisualization(self, vectorField, vectorFieldPos=None, f_vertices=None, f_faces=None, edgeColorFlag=0):
        if vectorField.ndim == 1:
            vectorField = vectorField.reshape(-1, 3)
        if vectorFieldPos is None:
            vectorFieldPos = self.baryCentersCalc()
        p = self.visualizeMesh(f_vertices, f_faces, edgeColorFlag)
        ax = plt.gca()
        ax.quiver(vectorFieldPos[:, 0], vectorFieldPos[:, 1], vectorFieldPos[:, 2],
                  vectorField[:, 0], vectorField[:, 1], vectorField[:, 2], length=0.1, normalize=True)
        return p

    def vectorFieldVisualization2(self, vectorField1, vectorField2, vectorFieldPos=None, f_vertices=None, f_faces=None, edgeColorFlag=0):
        p = self.visualizeMesh(f_vertices, f_faces, edgeColorFlag)
        ax = plt.gca()
        ax.quiver(vectorFieldPos[:, 0], vectorFieldPos[:, 1], vectorFieldPos[:, 2],
                  vectorField1[:, 0], vectorField1[:, 1], vectorField1[:, 2], color='b', length=0.1, normalize=True)
        ax.quiver(vectorFieldPos[:, 0], vectorFieldPos[:, 1], vectorFieldPos[:, 2],
                  vectorField2[:, 0], vectorField2[:, 1], vectorField2[:, 2], color='r', length=0.1, normalize=True)
        return p

    def visualizeDistances(self, u, x0, nisolines=0, urange=None, cam=None):
        if urange is None:
            urange = [u.min(), u.max()]
        p = self.visualizeMesh(u, edgeColorFlag=1, figFlag=1)
        ax = plt.gca()
        ax.set_xlim(urange)
        ax.scatter(self.vertices[x0, 0], self.vertices[x0, 1], self.vertices[x0, 2], c='r', s=100)
        if nisolines > 0:
            ax.tricontour(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2], self.faces, u, nisolines, cmap='jet')
        else:
            ax.set_cmap('jet')
        if cam is not None:
            ax.view_init(elev=cam[0], azim=cam[1])
        else:
            ax.view_init(elev=90, azim=0)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.tight_layout()
        return p

    def godf(self, n):
        M = self
        inner = M.inner_edges
        t1 = abs(M.e2t[inner, 0])
        t2 = abs(M.e2t[inner, 1])
        oe = -np.ones(M.nie)
        ze = np.zeros(M.nie)
        EV = M.vertices[M.edges[inner, 1]] - M.vertices[M.edges[inner, 0]]
        EV = self.normalize_vf(EV)
        IN1 = np.arctan2(np.sum(EV * M.F2[t1], axis=1), np.sum(EV * M.F1[t1], axis=1))
        IN2 = np.arctan2(np.sum(EV * M.F2[t2], axis=1), np.sum(EV * M.F1[t2], axis=1))
        PT = n * (IN2 - IN1)
        II = np.repeat(np.arange(M.nie), 8)
        JJ = np.concatenate((t1, t1 + M.nf, t1, t1 + M.nf, t2, t2 + M.nf, t2, t2 + M.nf))
        SS = np.concatenate((np.cos(PT), -np.sin(PT), np.sin(PT), np.cos(PT), oe, ze, ze, oe))
        CovD = csr_matrix((SS, (II, JJ)), shape=(2 * M.nie, 2 * M.nf))
        Ws = diags(np.repeat(np.sqrt(M.ea[inner]), 2))
        oph = Ws @ CovD
        op = oph.T @ oph
        return op, oph

    @staticmethod
    def normv(vf):
        return np.sqrt(np.sum(vf ** 2, axis=1))

    @staticmethod
    def normalize_vf(vf):
        nv = MeshClass.normv(vf)
        nnv = vf / nv[:, np.newaxis]
        nnv[nv < 1e-15] = 0
        return nnv

    @staticmethod
    def get_camera(ca=None):
        if ca is None:
            ca = plt.gca()
        cam = ca.view_init()
        return cam

    @staticmethod
    def set_camera(ca, cam):
        ca.view_init(elev=cam[0], azim=cam[1])

    @staticmethod
    def cotLaplacian(mesh, L23=None, L13=None, L12=None):
        X = mesh.vertices
        T = mesh.faces
        nv = X.shape[0]

        inputL = sum([L23 is not None, L13 is not None, L12 is not None])
        if inputL < 3:
            L1 = MeshClass.normv(X[T[:, 1]] - X[T[:, 2]])
            L2 = MeshClass.normv(X[T[:, 0]] - X[T[:, 2]])
            L3 = MeshClass.normv(X[T[:, 0]] - X[T[:, 1]])
        else:
            L1 = L23
            L2 = L13
            L3 = L12

        A1 = (L2**2 + L3**2 - L1**2) / (2 * L2 * L3)
        A2 = (L1**2 + L3**2 - L2**2) / (2 * L1 * L3)
        A3 = (L1**2 + L2**2 - L3**2) / (2 * L1 * L2)
        A = np.arccos(np.column_stack((A1, A2, A3)))

        I = np.concatenate((T[:, 0], T[:, 1], T[:, 2]))
        J = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))
        S = 0.5 * np.concatenate((1 / np.tan(A[:, 2]), 1 / np.tan(A[:, 0]), 1 / np.tan(A[:, 1])))
        In = np.concatenate((I, J, I, J))
        Jn = np.concatenate((J, I, I, J))
        Sn = np.concatenate((-S, -S, S, S))
        W = csr_matrix((Sn, (In, Jn)), shape=(nv, nv))

        if inputL < 3:
            M = MeshClass.mass_matrix_barycentric(mesh)
            A = M.sum(axis=1)
        else:
            M = MeshClass.mass_matrix_barycentric(mesh, L1, L2, L3)
            A = M.sum(axis=1)

        return W, A

    @staticmethod
    def mass_matrix_barycentric(mesh, L1=None, L2=None, L3=None):
        T = mesh.faces
        inputL = sum([L1 is not None, L2 is not None, L3 is not None])
        if inputL < 3:
            Ar = mesh.ta
        else:
            s = (L1 + L2 + L3) / 2
            Ar = np.sqrt(s * (s - L1) * (s - L2) * (s - L3))

        nv = mesh.nv

        I = np.concatenate((T[:, 0], T[:, 1], T[:, 2]))
        J = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))
        Mij = np.repeat(Ar / 12, 3)
        Mji = Mij
        Mii = np.repeat(Ar / 6, 3)
        In = np.concatenate((I, J, I))
        Jn = np.concatenate((J, I, I))
        Mn = np.concatenate((Mij, Mji, Mii))
        M = csr_matrix((Mn, (In, Jn)), shape=(nv, nv))

        return M
