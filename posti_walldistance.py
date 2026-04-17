#!/usr/bin/env python3

import argparse
import importlib

import h5py
import numpy as np

from posti.core import read_mesh


def walldistance(meshfile, geometry, N, lobatto, args):
    m = importlib.import_module(f"posti.walldistance.{geometry}")
    x = read_mesh(meshfile, N, lobatto)
    nElems = x.shape[0]
    d = np.zeros((nElems, N + 1, N + 1, N + 1))
    for iElem in range(nElems):
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    d[iElem, i, j, k] = m.walldistance(x[iElem, i, j, k], args)
        print(f"{iElem}/{nElems}", end="\r")
    print("")
    if meshfile.endswith("_mesh.h5"):
        wallfile = f"{meshfile[:-len('_mesh.h5')]}_walldistance.h5"
    else:
        wallfile = "walldistance.h5"
    f = h5py.File(wallfile, "w")
    f["walldistance"] = d
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile")
    parser.add_argument("geometry")
    parser.add_argument("N", type=int)
    parser.add_argument("--lobatto", action="store_true")
    parser.add_argument("--args", action="append", type=str)
    args = parser.parse_args()
    if args.args is None:
        args_dict = {}
    else:
        args_dict = dict(s.split("=", 1) for s in args.args)
    walldistance(args.meshfile, args.geometry, args.N, args.lobatto, args_dict)
