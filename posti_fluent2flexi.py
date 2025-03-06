#!/usr/bin/env python3

import argparse

import h5py
import numpy as np
import pandas as pd

from posti.core import read_mesh_state
from posti.swapmesh import SwapMesh


def fluent_to_flexi(fluentfile, meshfile, statefile):
    df = pd.read_csv(fluentfile, engine="python", sep=r",\s*")
    x_old = np.array(df[["x-coordinate", "y-coordinate"]])
    UPrim_old = np.array(
        df[
            [
                "density",
                "x-velocity",
                "y-velocity",
                "pressure",
                "turb-kinetic-energy",
                "specific-diss-rate",
            ]
        ]
    )
    swap = SwapMesh(x_old, UPrim_old)

    x, U, info = read_mesh_state(meshfile, statefile)
    nElems = info["nElems"]
    N = info["N"]

    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    rho, u, v, p, tke, omg = swap.query(x[iElem, i, j, k])
                    rhoE = p / (1.4 - 1.0) + 0.5 * rho * (u**2 + v**2) + rho * tke
                    g = np.sqrt(1.0 / (0.09 * max(omg, 1.0e-16)))
                    U[iElem, i, j, k] = np.array(
                        [rho, rho * u, rho * v, 0, rhoE, rho * tke, rho * g]
                    )
        print(f"{iElem}/{nElems}", end="\r")
    print(f"Max distance: {swap.max_distance()}")

    f = h5py.File(statefile, "r+")
    f["DG_Solution"][:] = U
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fluentfile")
    parser.add_argument("meshfile")
    parser.add_argument("statefile")
    args = parser.parse_args()
    fluent_to_flexi(args.fluentfile, args.meshfile, args.statefile)
