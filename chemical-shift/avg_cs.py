#!/usr/bin/python3
# This script is used to calculate the average of chemical shifts
# for each residues

import os, re, sys

import numpy as np

calc_frame = 500  # 构象个数！
mode = 's'

print(' Averaging chemical shifts ')

chemical = [[]]  # chemical[i][j][k], i for trajectory, j for frame, and k for atom

res_atom = []  # res_atom[i], i for name of residues and atom [A1_CA, A1_CB]

for i in range(calc_frame):

    chemical[0].append([])
    # sparta+
    tab_file = 'tmp_ab40_vae/drkN_%d.tab' % (i * 79 + 1)
    # if i == 1830:
    # 	continue
    # Check.FileExists(tab_file)
    tab = open(tab_file, 'r')

    for line in tab.readlines():
        if re.search('^[A-Z]', line) or re.search('^$', line):
            continue

        temp = re.split('\s+', line.strip('\n'))
        if i == 0:
            name = '%s%s_%s' % (temp[2], temp[1], temp[3])
            res_atom.append(name)

        chemical[0][i].append(float(temp[4]))

chemical = np.array(chemical)

chemical_avg = []
chemical_std = []
traj = 1
# print(calc_frame)

for k in range(len(res_atom)):
    avg = []
    if mode == 'm':  # the averaged value is the average of multiple trajectories
        for i in range(traj):
            temp_avg = 0
            for j in range(calc_frame):
                temp_avg += chemical[i][j][k]
            temp_avg /= calc_frame
            avg.append(temp_avg)
        chemical_avg.append(np.mean(avg))
        chemical_std.append(np.std(avg))

    elif mode == 's':  # the averaged value is the average of multiple blocks [4 blocks]
        blocks = 5
        for i in range(1):  # traj == 1 here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            temp_avg = 0
            count = 0
            for j in range(calc_frame):
                temp_avg += chemical[i][j][k]
                count += 1
                if count == calc_frame / blocks:
                    temp_avg /= count
                    avg.append(temp_avg)
                    temp_avg = 0
                    count = 0
        # if k == 4:
        # print(avg)
        chemical_avg.append(np.mean(avg))
        chemical_std.append(np.std(avg))
# print('%s %.4f %.4f' % (res_atom[k], chemical_avg[k], chemical_std[k])

# export result for each kind of atom
backbone = ['CA', 'CB', 'CO', 'N', 'HA', 'HN']

out_file = './RS1_10000_vae.dat'

out = open(out_file, 'w')

for k in range(len(res_atom)):
    res = res_atom[k].split('_')[0]
    atom = res_atom[k].split('_')[1]
    if atom == 'C':
        atom_name = 'CO'
    elif atom == 'H':
        atom_name = 'HN'
    elif atom == 'HA2':
        atom_name = 'HA'
    else:
        atom_name = atom
    if atom_name in backbone:
        out.write('%s\t%s\t%.4f\t%.4f\n' % (res, atom_name, chemical_avg[k], chemical_std[k]))

out.close()
