import os
from math import dist
import matplotlib.pyplot as plt
import numpy as np

r_alpha_length = 0
r_310_length = 0
sheet_length = 0
total_length = 0
sheet_length_21 = 0

total_residues = {}
r_alpha_residues = {}
sheet_residues = {}

total_distance_list = []
r_alpha_distance_list=[]


for file in os.listdir("Supplementary"):
    print(r_alpha_residues)
    residueID_residueName = {}
    N_residueID_coords = {}
    O_residueID_coords = {}
    with open(f"Supplementary/{file}") as f:
        for line in f.readlines():
            if line.startswith("ATOM"):
                chain_id = line[21]
                residue_id = int(line[22:26])
                aa_name = line[17:20].strip()
                if chain_id not in residueID_residueName.keys():
                    residueID_residueName[chain_id] = {}
                    N_residueID_coords[chain_id] = {}
                    O_residueID_coords[chain_id] = {}

                if aa_name == 3: # remove some entry like 3q0n binding with rna
                    residueID_residueName[chain_id][residue_id] = aa_name
                    if line[12:16].strip() == 'N':
                        N_residueID_coords[chain_id][residue_id] = [float(line[30:38]), float(line[38:46]),
                                                                    float(line[46:54])]
                    if line[12:16].strip() == 'O':
                        O_residueID_coords[chain_id][residue_id] = [float(line[30:38]), float(line[38:46]),
                                                                    float(line[46:54])]

    for v in residueID_residueName.values():
        for aa in v.values():
            if aa not in total_residues.keys():
                total_residues[aa] = 0
            total_residues[aa] += 1

    with open(f"Supplementary/{file}") as f:
        for line in f.readlines():
            if line.startswith("HELIX"):
                chainID = line[19]
                helix_class = int(line[38:40])
                if helix_class == 1:  # right alpha
                    r_alpha_length += int(line[71:76])  # sum the length od helix
                    for i in range(int(line[21:25]), int(line[33:37]) + 1):

                        if i in residueID_residueName[chainID].keys():
                            if residueID_residueName[chainID][i] not in r_alpha_residues.keys():
                                r_alpha_residues[residueID_residueName[chainID][i]] = 0
                            r_alpha_residues[residueID_residueName[chainID][i]] += 1

                            # distance calculating
                            if i + 4 in residueID_residueName[chainID].keys() and i + 4 <= int(line[33:37]):
                                distance = dist(N_residueID_coords[chainID][i + 4],
                                                O_residueID_coords[chainID][i])
                                r_alpha_distance_list.append(distance)

                if helix_class == 5:  # 310 helix
                    r_310_length += int(line[71:76])

            if line.startswith("SHEET"):
                sheet_length_21 += int(line[33:37]) - int(line[22:26]) +1
                chainID=line[21]
                for i in range(int(line[22:26]), int(line[33:37]) + 1):
                    if i in residueID_residueName[chainID].keys():
                        if residueID_residueName[chainID][i] not in sheet_residues.keys():
                            sheet_residues[residueID_residueName[chainID][i]] = 0
                        sheet_residues[residueID_residueName[chainID][i]] += 1

    for chainID, value in O_residueID_coords.items():
        for residueID in value.keys():
            if residueID + 4 in value.keys():
                distance = dist(N_residueID_coords[chainID][residueID + 4], O_residueID_coords[chainID][residueID])
                total_distance_list.append(distance)

sheet_length = sum([v for v in sheet_residues.values()])
total_length += sum([v for v in total_residues.values()])
print(r_alpha_length / total_length)
print(r_310_length / total_length)
print(sheet_length / total_length)
print(total_length)
print('sum of length right alpha', r_alpha_length)
print('matched residues right alpha', sum([v for v in r_alpha_residues.values()]))
print('matched residues sheet', sheet_length)
print('sheet by position', sheet_length_21)

for k in total_residues.keys():
    print(k, "\t", r_alpha_residues[k], "\t", r_alpha_residues[k] / total_residues[k])
print()
for k in total_residues.keys():
    print(k, "\t", sheet_residues[k], "\t", sheet_residues[k] / total_residues[k])

bins = np.linspace(min(total_distance_list), max(total_distance_list), 100)
plt.hist(total_distance_list, bins=bins)
plt.xlabel('Distance')
plt.ylabel('Frequency')
plt.title('all structures distance')
plt.savefig("hist_all_distance.png")
plt.clf()
plt.hist(r_alpha_distance_list, bins=bins)
plt.xlabel('Distance')
plt.ylabel('Frequency')
plt.title('right alpha helix distance')
plt.savefig("hist_r_alpha_helix_distance.png")