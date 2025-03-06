import numpy as np
import scipy.special
import h5py


def vandermonde(input, output):
    m = len(output)
    n = len(input)
    A = np.ones((m, n))
    for i in range(m):
        for j in range(n):
            for l in range(n):
                if l != j:
                    A[i, j] *= (output[i] - input[l]) / (input[j] - input[l])
    return A


def Ngeo_to_N_vdm(Ngeo, N):
    x0 = np.array([-1 + i * 2 / Ngeo for i in range(Ngeo + 1)])
    x1 = np.array([np.cos(i / Ngeo * np.pi) for i in range(Ngeo + 1)])
    x2 = np.array([np.cos(i / N * np.pi) for i in range(N + 1)])
    x3, _ = scipy.special.roots_legendre(N + 1)
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


def build_coords(NodeCoords, nElems, Ngeo, N):
    x = np.zeros((nElems, N + 1, N + 1, N + 1, 3))
    A = Ngeo_to_N_vdm(Ngeo, N)
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


def read_mesh(meshfile, N):
    f = h5py.File(meshfile)
    Ngeo = int(f.attrs["Ngeo"][0])
    nElems = int(f.attrs["nElems"][0])
    NodeCoords = f["NodeCoords"]
    x = build_coords(NodeCoords, nElems, Ngeo, N)
    f.close()
    return x


def read_mesh_state(meshfile, statefile, upsampling=0):
    s = h5py.File(statefile)
    nElems = s["DG_Solution"].shape[0]
    N = s["DG_Solution"].shape[1] - 1
    nVars = s["DG_Solution"].shape[4]
    U = np.zeros(s["DG_Solution"].shape)
    U[:] = s["DG_Solution"]
    s.close()
    x = read_mesh(meshfile, N)
    if upsampling == 0:
        info = {"nElems": nElems, "N": N, "nVars": nVars}
        return x, U, info
    else:
        if upsampling < 0:
            N_new = int(-upsampling * (N + 1)) - 1
        else:
            N_new = upsampling
        assert N_new > N
        x_new = np.zeros((nElems, N_new + 1, N_new + 1, N_new + 1, 3))
        U_new = np.zeros((nElems, N_new + 1, N_new + 1, N_new + 1, nVars))
        x0, _ = scipy.special.roots_legendre(N + 1)
        x1, _ = scipy.special.roots_legendre(N_new + 1)
        A = vandermonde(x0, x1)
        for iElem in range(nElems):
            x_new[iElem] = change_basis(x[iElem], A)
            U_new[iElem] = change_basis(U[iElem], A)
        info = {"nElems": nElems, "N": N_new, "nVars": nVars}
        return x_new, U_new, info
