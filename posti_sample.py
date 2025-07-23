#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd

from posti.core import read_mesh_state


def sample(meshfile, statefile, nx, ny, nz, x_sample, outputfile):
    x, U, info = read_mesh_state(meshfile, statefile, readmean=True)
    nElems = info["nElems"]
    N = info["N"]
    nVars = info["nVars"]
    varNames = info["varNames"]

    x = x.reshape((nElems * (N + 1) ** 3, 3))
    U = U.reshape((nElems * (N + 1) ** 3, nVars))
    data = {
        "x": x[:, 0],
        "y": x[:, 1],
        "z": x[:, 2],
    }
    for i, varName in enumerate(varNames):
        data[varName] = U[:, i]

    U2 = np.zeros((nx, ny, 3 + nVars))

    df = pd.DataFrame(data=data)
    df.sort_values(by="x", inplace=True)
    for i in range(nx):
        df_x = df.iloc[i * ny * nz : (i + 1) * ny * nz].copy()
        df_x.sort_values(by="y", inplace=True)
        for j in range(ny):
            df_y = df_x.iloc[j * nz : (j + 1) * nz].copy()
            UU = df_y.to_numpy()
            U2[i, j] = np.average(UU, axis=0)

    for i in range(nx):
        if U2[i, 0, 0] >= x_sample:
            kx = (U2[i, 0, 0] - x_sample) / (U2[i, 0, 0] - U2[i - 1, 0, 0])
            UU = kx * U2[i - 1, :, :] + (1 - kx) * U2[i, :, :]
            columns = ["x", "y", "z"]
            for varName in varNames:
                columns.append(varName.decode("utf-8"))
            output = pd.DataFrame(data=UU.reshape((ny, 3 + nVars)), columns=columns)
            output.drop(columns=["x", "z"], inplace=True)
            output.to_csv(outputfile, index=False)
            break

    # columns = ["x", "y", "z"]
    # for varName in varNames:
    #     columns.append(varName.decode("utf-8"))
    # output = pd.DataFrame(data=U2.reshape((nx * ny, 3 + nVars)), columns=columns)
    # output.drop(columns=["z"], inplace=True)
    # output.to_csv(outputfile, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile")
    parser.add_argument("statefile")
    parser.add_argument("nx", type=int)
    parser.add_argument("ny", type=int)
    parser.add_argument("nz", type=int)
    parser.add_argument("x", type=float)
    parser.add_argument("outputfile")
    args = parser.parse_args()
    sample(
        args.meshfile,
        args.statefile,
        args.nx,
        args.ny,
        args.nz,
        args.x,
        args.outputfile,
    )
