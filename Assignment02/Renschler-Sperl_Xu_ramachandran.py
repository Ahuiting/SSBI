# Author: Huiting Xu and Desirée Renschler-Sperl
import argparse
import os

from Bio.PDB import *
import numpy as np
import numpy.linalg as la
import math
import matplotlib.pyplot as plt


def create_parser():
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', nargs='+', help='PATH_TO_PDB_FILE')
    p.add_argument('-o', '--output', help='PATH_FOR_PDF_OUTPUT')
    return p.parse_args()


def get_N_CA_C_list(pdbfile):
    parser = PDBParser()
    structure = parser.get_structure(pdbfile[:4], pdbfile)
    atom_dict = {}
    for chain in structure.get_chains():
        atom_list = []
        atom_dict[chain.get_id()] = atom_list
        for residue in chain:
            atom_residue_list = []  # N, CA, C
            for atom in residue:
                if atom.get_name() == 'N' or atom.get_name() == 'CA' or atom.get_name() == 'C':
                    atom_residue_list.append(atom)
            if len(atom_residue_list) != 0:
                atom_list.append(atom_residue_list)
    return atom_dict


def get_angle(v1, v2):
    cos = np.dot(v1, v2)
    sin = la.norm(np.cross(v1, v2))
    return 180 / math.pi * np.arctan2(sin, cos)


def get_normal_vector(point1, point2, point3):
    v_12 = point2 - point1
    v_13 = point3 - point1
    n_vector = np.cross(v_12, v_13)
    return n_vector / np.linalg.norm(n_vector)


def get_phi_psi_list(atom_dict):
    phi_list = []
    psi_list = []
    for v in atom_dict.values():
        for i in range(1, len(v) - 1):
            prevC = v[i - 1][2].get_coord()
            currN = v[i][0].get_coord()
            currCa = v[i][1].get_coord()
            currC = v[i][2].get_coord()
            nextN = v[i + 1][0].get_coord()

            l_plane = get_normal_vector(prevC, currN, currCa)
            m_plane = get_normal_vector(currN, currCa, currC)
            r_plane = get_normal_vector(currCa, currC, nextN)

            phi = get_angle(l_plane, m_plane)
            psi = get_angle(m_plane, r_plane)

            if np.dot(prevC - currCa, m_plane) > 0:
                phi = -phi
            if np.dot(nextN - currCa, m_plane) < 0:
                psi = -psi

            phi_list.append(phi)
            psi_list.append(psi)
    return phi_list, psi_list


if __name__ == '__main__':
    'try -i 5ire.pdb 1igt.pdb -o RamachandranMaps'
    args = create_parser()
    file_list = args.input
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    for file in file_list:
        phi_list, psi_list = get_phi_psi_list(get_N_CA_C_list(file))
        # plot
        plt.clf()
        plt.scatter(phi_list, psi_list, s=1)
        plt.xlabel('phi')
        plt.ylabel('psi')
        plt.title(f'{file[:4]} Ramachandran Map')
        plt.axhline(y=0, color='k', linestyle='-', linewidth=1)
        plt.axvline(x=0, color='k', linestyle='-', linewidth=1)
        plt.savefig(f'{args.output}/{file[:4]}_Ramachandran_Map')
