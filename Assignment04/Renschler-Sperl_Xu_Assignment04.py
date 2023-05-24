# Author: Huiting Xu and Desirée Renschler-Sperl
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
import math
from scipy.constants import epsilon_0, elementary_charge
import seaborn as sns
import matplotlib.pyplot as plt


def get_residue_list(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure(pdb_file[:4], pdb_file)

    residue_list = []
    for residue in structure[0].get_residues():
        residue_id = residue.get_id()[1]
        c_coords, o_coords, n_coords, h_coords = None, None, None, None
        for atom in residue:
            if atom.get_name() == "C":
                c_coords = atom.get_coord()
            elif atom.get_name() == "O":
                o_coords = atom.get_coord()
            elif atom.get_name() == "N":
                n_coords = atom.get_coord()
            elif atom.get_name() == "H":
                h_coords = atom.get_coord()
        residue_list.append((residue_id, c_coords, o_coords, n_coords, h_coords))
    return residue_list


def calculate_energy(c_coords, o_coords, n_coords, h_coords):
    q1 = 0.42 * elementary_charge
    q2 = 0.20 * elementary_charge

    d_ON = math.dist(o_coords, n_coords)
    d_CH = math.dist(c_coords, h_coords)
    d_OH = math.dist(o_coords, h_coords)
    d_CN = math.dist(c_coords, n_coords)

    energy = (q1 * q2 / (4 * math.pi * epsilon_0)) * (1 / d_ON + 1 / d_CH - 1 / d_OH - 1 / d_CN)
    return energy


if __name__ == '__main__':
    residue_list=get_residue_list("5jxv.pdb")

    matrix_size = len(residue_list)
    dssp_matrix = np.zeros([matrix_size, matrix_size])
    for i in range(matrix_size):
        for j in range(matrix_size):
            if i != j:
                _, c_coords, o_coords, _, _ = residue_list[i]
                _, _, _, n_coords, h_coords = residue_list[j]
                energy = calculate_energy(c_coords, o_coords, n_coords, h_coords)
                dssp_matrix[i][j] = energy

    df = pd.DataFrame(dssp_matrix)
    df.index = [i[0] for i in residue_list]
    df.to_csv("dssp_matrix.tsv", sep="\t", header=[i[0] for i in residue_list], index=True)

    sns.heatmap(dssp_matrix, vmax=sum(dssp_matrix.flatten() / 56), vmin=sorted(dssp_matrix.flatten())[-1])
    plt.title("DSSP Energy Matrix")
    plt.xlabel("Residue Index j")
    plt.ylabel("Residue Index i")
    plt.savefig("Heatmap.png")
