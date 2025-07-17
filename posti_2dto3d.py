#!/usr/bin/env python3

import argparse

import h5py

from posti.core import read_mesh_state
from posti.swapmesh import SwapMesh


def from2dto3d(meshfile_old, statefile_old, meshfile, statefile):
    x_old, U_old, info_old = read_mesh_state(meshfile_old, statefile_old)
    nElems_old = info_old["nElems"]
    N_old = info_old["N"]
    nVars_old = info_old["nVars"]

    swap = SwapMesh(
        x_old.reshape((nElems_old * (N_old + 1) ** 3, 3))[:, 0:2],
        U_old.reshape((nElems_old * (N_old + 1) ** 3, nVars_old)),
    )

    x, U, info = read_mesh_state(meshfile, statefile)
    nElems = info["nElems"]
    N = info["N"]
    nVars = info["nVars"]

    # assert N_old == N
    # assert nVars_old == nVars

    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    U[iElem, i, j, k] = swap.query(x[iElem, i, j, k], from2dto3d=True, N = N + 1)
        print(f"{iElem}/{nElems}", end="\r")
    print("")
    print(f"Max distance: {swap.max_distance()}")

    f = h5py.File(statefile, "r+")
    f["DG_Solution"][:] = U
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile_old")
    parser.add_argument("statefile_old")
    parser.add_argument("meshfile")
    parser.add_argument("statefile")
    args = parser.parse_args()
    from2dto3d(args.meshfile_old, args.statefile_old, args.meshfile, args.statefile)
