import argparse
import numpy as np
from Bio import PDB


def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', nargs=2, help='PDB_FILE')
    return p.parse_args()


def calculate_rmsd(file1, file2):
    parser = PDB.PDBParser()
    file1_structure = parser.get_structure('file1', file1)
    file2_structure = parser.get_structure('file2', file2)

    file1_atoms = np.array([i.get_coord() for i in file1_structure.get_atoms()])
    file2_atoms = np.array([i.get_coord() for i in file2_structure.get_atoms()])

    rmsd = np.sqrt(np.mean(np.sum((file1_atoms - file2_atoms) ** 2, axis=1)))
    return rmsd


if __name__ == '__main__':
    args = create_parser()
    rmsd_value = calculate_rmsd(args.input[0], args.input[1])
    print("RMSD value:", rmsd_value)
