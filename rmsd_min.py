import sys
from copy import deepcopy

import biobox as bb
import numpy as np

infile = sys.argv[1]
encoding_dim = 2
BATCH_SIZE = 64
EPOCHS = 50
split = sys.argv[2]


def cal_rmsd(input_data, output, test_size):
    test = deepcopy(input_data)
    examine = deepcopy(output)
    re = []
    for j in range(0, test_size):
        test.add_xyz(examine[j])
    for i in range(0, test_size, 3):
        for j in range(0, test_size, 3):
            val = test.rmsd(i, test_size + j)
            # print("({}, {})".format(i, j))
            re.append(val)
    return re


# original conformations
M = bb.Molecule()
M.import_pdb("./data/%s" % infile.replace(".pdb", ".aligned.pdb"))
idx = M.atomselect("*", "*", ["CA", "CB", "C", "N", "CG", "CD", "NE2", "CH3", "O", "OE1"], get_index=True)[1]
crds = M.coordinates[:, idx]
indices = np.random.permutation(len(crds))

# structures of test set for later comparison
indicestest = np.sort(indices[int(split * M.coordinates.shape[0]), M.coordinates.shape[0], 1])
M2 = M.get_subset(idxs=idx, conformations=indicestest)

M1 = bb.Molecule()
M1.import_pdb("./data/%s" % infile.replace(".pdb", ".generated.pdb"))
idx = M1.atomselect("*", "*", ["CA", "CB", "C", "N", "CG", "CD", "NE2", "CH3", "O", "OE1"], get_index=True)[1]
M3 = M1.get_subset(idxs=idx)

comp_test = cal_rmsd(M2, M3, test_size=(M.coordinates.shape[0] - int(split * M.coordinates.shape[0])))
comp_test = np.array(comp_test)
print("> average rmsd: %.3f" % np.mean(comp_test))
np.savetxt("./data/rmsd.dat")
