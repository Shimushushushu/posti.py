#!/usr/bin/env python3

import argparse

import h5py
import numpy as np

from posti.core import read_mesh
from posti.walldistance import channel, faith, flatplate


def walldistance(meshfile, geometry, N):
    x = read_mesh(meshfile, N)
    nElems = x.shape[0]
    d = np.zeros((nElems, N + 1, N + 1, N + 1))
    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    if geometry == "channel":
                        d[iElem, i, j, k] = channel.walldistance(x[iElem, i, j, k])
                    elif geometry == "flatplate":
                        d[iElem, i, j, k] = flatplate.walldistance(x[iElem, i, j, k])
                    elif geometry == "faith":
                        d[iElem, i, j, k] = faith.walldistance(x[iElem, i, j, k])
        print(f"{iElem}/{nElems}", end="\r")
    print("")
    f = h5py.File("walldistance.h5", "w")
    f["walldistance"] = d
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile")
    parser.add_argument("geometry")
    parser.add_argument("N", type=int)
    args = parser.parse_args()
    walldistance(args.meshfile, args.geometry, args.N)
