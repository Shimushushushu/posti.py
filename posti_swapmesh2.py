#!/usr/bin/env python3

import argparse

import h5py

from posti.core import read_mesh_state
from posti.swapmesh import SwapMesh


def swapmesh2(meshfile_old, statefile_old, meshfile, statefile, upsampling):
    x_old, U_old, info_old = read_mesh_state(meshfile_old, statefile_old, upsampling)
    nElems_old = info_old["nElems"]
    N_old = info_old["N"]
    nVars = info_old["nVars"]
    x_old.resize((nElems_old * (N_old + 1) ** 3, 3))
    U_old.resize((nElems_old * (N_old + 1) ** 3, nVars))
    swap = SwapMesh(x_old, U_old)

    x, U, info = read_mesh_state(meshfile, statefile)
    nElems = info["nElems"]
    N = info["N"]
    assert nVars == info["nVars"]

    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    U[iElem, i, j, k] = swap.query(x[iElem, i, j, k])
        print(f"{iElem}/{nElems}", end="\r")
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
    parser.add_argument(
        "-N",
        "--upsampling",
        type=float,
        default=-1.5,
        help="Upsample the old state file.  If N = 0, no upsampling.  If N > 0, upsample to the given N.  If N < 0, upsample to (-N * N_old), where N_old is the original N.",
    )
    args = parser.parse_args()
    swapmesh2(
        args.meshfile_old,
        args.statefile_old,
        args.meshfile,
        args.statefile,
        args.upsampling,
    )
