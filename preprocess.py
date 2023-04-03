import sys
import time
import os
from copy import deepcopy

import biobox as bb
import numpy as np


def main():
    x_train_orig, x_test_orig = load_pdb(sys.argv[1], split=float(sys.argv[2]), aligned=False)

    if not os.path.exists("./data/"):
        os.makedirs("./data/")
    np.savetxt("./data/%s" % (sys.argv[1].replace(".pdb", ".train.dat")), x_train_orig)
    np.savetxt("./data/%s" % (sys.argv[1].replace(".pdb", "_test.dat")), x_test_orig)


def load_pdb(filename, split=0.2, aligned=False):
    print("\n> loading molecule ", filename)
    M = bb.Molecule()
    M.import_pdb(filename)

    # atom selection for reconstruction
    idx = M.atomselect("*", "*", ["CA", "CB", "C", "N", "CG", "CD", "NE2", "CH3", "O", "OE1"], get_index=True)[1]
    crds = M.coordinates[:, idx]
    print("> molecule has %s conformations" % M.coordinates.shape[0])
    print("> selected %s atoms" % len(idx))

    # align all structures, and center to origin
    if not aligned:
        print("> align conformations ...")
        time0 = time.time()
        align(M)
        print("  alignment completed, takes %.2f s" % (time.time()-time0))

    # save structures of train and test set for later comparison
    indicestrain = np.arange(0, int(split * M.coordinates.shape[0]), 1)
    indicestest = np.arange(int(split * M.coordinates.shape[0]), M.coordinates.shape[0], 1)

    x_train = deepcopy(crds[indicestrain])
    x_train = x_train.reshape(x_train.shape[0], np.prod(x_train.shape[1:]))

    x_test = deepcopy(crds[indicestest])
    x_test = x_test.reshape(x_test.shape[0], np.prod(x_test.shape[1:]))

    return x_train, x_test


def align(M):
    idx1 = M.atomselect("*", "*", "CA", get_index=True)[1]
    M.rmsd_one_vs_all(0, points_index=idx1, align=True)

    for pos in range(len(M.coordinates)):
        M.set_current(pos)
        M.center_to_origin()
    M.set_current(0)
    M.write_pdb("./data/%s" % sys.argv[1].replace(".pdb", ".aligned.pdb"))


if __name__ == "__main__":
    main()
