# Author: Huiting Xu and Desirée Renschler-Sperl
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

    file1_name = file1[:-4]
    file2_name = file2[:-4]

    cmd.align(file1_name, file2_name, object='align')

    file1_atoms_coords = np.asarray(cmd.get_model(f'align and {file1_name}').get_coord_list())
    file2_atoms_coords = np.asarray(cmd.get_model(f'align and {file2_name}').get_coord_list())

    rmsd = np.sqrt(np.mean(np.sum((file1_atoms_coords - file2_atoms_coords) ** 2, axis=1)))
    return rmsd


if __name__ == '__main__':
    args = create_parser()
    rmsd_value = calculate_rmsd(args.input[0], args.input[1])
    print("RMSD value:", rmsd_value)
