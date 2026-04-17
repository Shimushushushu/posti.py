import h5py
import numpy as np
from scipy.interpolate import lagrange
from scipy.special import legendre, roots_legendre


def vandermonde(input, output):
    m = len(output)
    n = len(input)
    A = np.ones((m, n))
    for j in range(n):
        f = np.zeros(n)
        f[j] = 1
        l = lagrange(input, f)
        A[:, j] = l(output)
    return A


def roots_legendre_lobatto(n):
    if n == 2:
        x = np.array([-1.0, 1.0])
        w = np.array([1.0, 1.0])
        return x, w
    else:
        Pn = legendre(n - 1)
        dPn = Pn.deriv()
        x = np.sort(np.append(dPn.roots, [-1.0, 1.0]))
        w = 2.0 / (n * (n - 1) * (Pn(x) ** 2))
        return x, w


def nodal_points(N, lobatto=False):
    if lobatto:
        return roots_legendre_lobatto(N + 1)
    else:
        return roots_legendre(N + 1)


def Ngeo_to_N_vdm(Ngeo, N, lobatto=False):
    x0 = np.array([-1 + i * 2 / Ngeo for i in range(Ngeo + 1)])
    x1 = np.array([np.cos(i / Ngeo * np.pi) for i in range(Ngeo + 1)])
    x2 = np.array([np.cos(i / N * np.pi) for i in range(N + 1)])
    x3, _ = nodal_points(N, lobatto)
    A1 = vandermonde(x0, x1)
    A2 = vandermonde(x1, x2)
    A3 = vandermonde(x2, x3)
    return A3 @ A2 @ A1


def change_basis(U, A):
    N1, N2, N3, Nv = U.shape
    No, Ni = A.shape
    assert N1 == N2 == N3 == Ni
    U1 = np.zeros((Ni, Ni, No, Nv))
    U2 = np.zeros((Ni, No, No, Nv))
    U3 = np.zeros((No, No, No, Nv))
    for i in range(Ni):
        for j in range(Ni):
            for k in range(No):
                for l in range(Ni):
                    U1[i, j, k] += A[k, l] * U[i, j, l]
    for i in range(Ni):
        for j in range(No):
            for l in range(Ni):
                for k in range(No):
                    U2[i, j, k] += A[j, l] * U1[i, l, k]
    for i in range(No):
        for l in range(Ni):
            for j in range(No):
                for k in range(No):
                    U3[i, j, k] += A[i, l] * U2[l, j, k]
    return U3


def build_coords(NodeCoords, nElems, Ngeo, N, lobatto=False):
    x = np.zeros((nElems, N + 1, N + 1, N + 1, 3))
    A = Ngeo_to_N_vdm(Ngeo, N, lobatto)
    iNode = 0
    for iElem in range(nElems):
        x_geo = np.zeros((Ngeo + 1, Ngeo + 1, Ngeo + 1, 3))
        for i in range(Ngeo + 1):
            for j in range(Ngeo + 1):
                for k in range(Ngeo + 1):
                    x_geo[i, j, k] = NodeCoords[iNode]
                    iNode += 1
        x[iElem] = change_basis(x_geo, A)
    return x


def read_mesh(meshfile, N, lobatto=False):
    f = h5py.File(meshfile)
    Ngeo = int(f.attrs["Ngeo"][0])
    nElems = int(f.attrs["nElems"][0])
    NodeCoords = f["NodeCoords"]
    x = build_coords(NodeCoords, nElems, Ngeo, N, lobatto)
    f.close()
    return x


def read_mesh_state(meshfile, statefile):
    s = h5py.File(statefile)
    lobatto = s.attrs["NodeType"][0].decode() == "GAUSS-LOBATTO"
    is2d = s.attrs["Dimension"][0] == 2
    file_type = s.attrs["File_Type"][0].decode()
    if file_type == "State":
        dataset_name = "DG_Solution"
        varNames = list(map(lambda x: x.decode(), s.attrs["VarNames"]))
    elif file_type == "TimeAvg":
        dataset_name = "Mean"
        varNames = list(map(lambda x: x.decode(), s.attrs["VarNames_Mean"]))
    else:
        print(f"error: not known file type: {file_type}")
        return None
    nElems = s[dataset_name].shape[0]
    N = s[dataset_name].shape[3] - 1
    nVars = s[dataset_name].shape[4]
    U = np.zeros((nElems, N + 1, N + 1, N + 1, nVars))
    if is2d:
        for i in range(N + 1):
            U[:, i : i + 1, :, :, :] = s[dataset_name]
    else:
        U[:] = s[dataset_name]
    s.close()
    x = read_mesh(meshfile, N, lobatto)
    info = {
        "nElems": nElems,
        "N": N,
        "nVars": nVars,
        "varNames": varNames,
        "is2d": is2d,
        "islobatto": lobatto,
    }
    return x, U, info
