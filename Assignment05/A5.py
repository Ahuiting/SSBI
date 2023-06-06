import argparse
import cmd

import numpy as np
from pymol import cmd


def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', nargs=2, help='PDB_FILE')
    return p.parse_args()


def calculate_rmsd(file1, file2):

    cmd.load(file1)
    cmd.load(file2)

    cmd.align('5uo8', '6pp4', object='align')
    alignment = cmd.get_raw_alignment("align")

    atom_dist = []
    for idx1, idx2 in alignment:
        atom1_coord = np.array(cmd.get_coords(f'{idx1[0]} and id {idx1[1]}'), dtype=float).flatten()
        atom2_coord = np.array(cmd.get_coords(f'{idx2[0]} and id {idx2[1]}'), dtype=float).flatten()
        if len(atom1_coord) == 3 and len(atom2_coord) == 3:
            dist = sum((atom1_coord[i] - atom2_coord[i]) ** 2 for i in range(len(atom1_coord)))
            atom_dist.append(dist)

    rmsd = np.sqrt(np.mean(atom_dist))
    return rmsd


if __name__ == '__main__':
    args = create_parser()
    rmsd_value = calculate_rmsd(args.input[0], args.input[1])
    print("RMSD value:", rmsd_value)
