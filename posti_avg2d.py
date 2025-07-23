#!/usr/bin/env python3

import argparse
import os

import h5py
import numpy as np

from posti.core import read_mesh_state
from posti.swapmesh import SwapMesh


def avg2d(meshfile, statefile, nDOFsSpanwise):
    os.system(f"cp {statefile} {statefile}.bak")
    x, U, info = read_mesh_state(meshfile, statefile, readmean=True)
    nElems = info["nElems"]
    N = info["N"]
    nVars = info["nVars"]

    swap = SwapMesh(
        x.reshape((nElems * (N + 1) ** 3, 3))[:, 0:2],
        U.reshape((nElems * (N + 1) ** 3, nVars)),
    )

    Uavg = np.zeros_like(U)
    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    Uavg[iElem, i, j, k] = swap.query(
                        x[iElem, i, j, k], from2dto3d=True, N=nDOFsSpanwise
                    )
        print(f"{iElem}/{nElems}", end="\r")
    print("")
    print(f"Max distance: {swap.max_distance()}")

    f = h5py.File(statefile, "r+")
    if "Mean" in f.keys():
        f["Mean"][:] = Uavg
    else:
        f["DG_Solution"][:] = Uavg
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile")
    parser.add_argument("statefile")
    parser.add_argument("nDOFsSpanwise", type=int)
    args = parser.parse_args()
    avg2d(args.meshfile, args.statefile, args.nDOFsSpanwise)
