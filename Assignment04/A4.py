import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
import math
from scipy.constants import epsilon_0, elementary_charge
import seaborn as sns
import matplotlib.pyplot as plt


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
    parser = PDBParser()
    structure = parser.get_structure("5JXV", "5jxv.pdb")

    residue_info = []
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
        residue_info.append((residue_id, c_coords, o_coords, n_coords, h_coords))

    matrix_size = len(residue_info)
    dssp_matrix = np.zeros([matrix_size, matrix_size])
    for i in range(matrix_size):
        for j in range(matrix_size):
            if i != j:
                _, c_coords, o_coords, _, _ = residue_info[i]
                _, _, _, n_coords, h_coords = residue_info[j]
                energy = calculate_energy(c_coords, o_coords, n_coords, h_coords) * 100000000000000000
                dssp_matrix[i][j] = energy
                dssp_matrix[j][i] = energy

    df = pd.DataFrame(dssp_matrix)
    df.index = [i[0] for i in residue_info]
    df.to_csv("dssp_matrix.tsv", sep="\t", header=[i[0] for i in residue_info], index=True)

    sns.heatmap(dssp_matrix,vmin=sum(dssp_matrix.flatten()/56),vmax=sorted(dssp_matrix.flatten())[-1])
    plt.title("DSSP Energy Matrix")
    plt.xlabel("Residue Index j")
    plt.ylabel("Residue Index i")
    plt.show()
