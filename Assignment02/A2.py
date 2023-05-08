# Author: Huiting Xu and Desirée Renschler-Sperl
from Bio.PDB import *
import numpy as np
import numpy.linalg as la
import math
import matplotlib.pyplot as plt

parser = PDBParser()
structure = parser.get_structure('1igt', "1igt.pdb")
atom_dict = {}
for model in structure:
    for chain in model:
        atom_list = []
        atom_dict[chain.get_id()] = atom_list
        for residue in chain:
            atom_residue_list = []  # 'N, CAC C
            for atom in residue:
                if atom.get_name() == 'N' or atom.get_name() == 'CA' or atom.get_name() == 'C':
                    atom_residue_list.append(atom)
            if len(atom_residue_list) != 0:
                atom_list.append(atom_residue_list)


def get_cross_product(v1, v2):
    new_v = [v1[1] * v2[2] - v1[2] * v2[1],
             v1[2] * v2[0] - v1[0] * v2[2],
             v1[0] * v2[1] - v1[1] * v2[0]]
    return new_v


def get_angle(v1, v2):
    cos = np.dot(v1, v2)
    sin = la.norm(np.cross(v1, v2))
    return  180/math.pi * np.arctan2(sin, cos)


phi_list=[]
psi_list=[]
for v in atom_dict.values():
    for i in range(1, len(v) - 1):
        prevC_currCa = v[i - 1][2].get_coord() - v[i][1].get_coord()
        currN_currCa = v[i][0].get_coord() - v[i][1].get_coord()
        currC_currCa = v[i][2].get_coord() - v[i][1].get_coord()
        nextN_currCa = v[i + 1][0].get_coord() - v[i][1].get_coord()

        l_plane = get_cross_product(prevC_currCa, currN_currCa)
        m_plane = get_cross_product(currC_currCa, currN_currCa)
        r_plane = get_cross_product(currC_currCa, nextN_currCa)

        phi=get_angle(l_plane, m_plane)
        psi = get_angle(m_plane, r_plane)

        if np.dot(prevC_currCa,m_plane)>0:
            phi=-phi
        if np.dot(nextN_currCa, m_plane) <0:
            psi = -psi
        phi_list.append(phi)
        psi_list.append(psi)



plt.scatter(phi_list,psi_list)
plt.show()




