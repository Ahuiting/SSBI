# Author: Huiting Xu and Desirée Renschler-Sperl
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from re import findall as refindall
from matplotlib import pyplot as plt

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input',  help='fasta_file')
    return p.parse_args()

def trypsin_digestion(sequence):
    peptide_list = []
    peptide = ""
    for index, aa in enumerate(sequence):
        peptide += aa
        if (aa == "K" or aa == "R") and (index + 1 < len(sequence)) and (sequence[index + 1] != "P"):
            peptide_list.append(peptide)
            peptide = ""
    # the last fragment
    if peptide:
        peptide_list.append(peptide)

    return peptide_list


# calculate mass for all the amino acids
def prepare_masses_dict():
    amino_acids = {
        "I": "C6H13NO2",
        "L": "C6H13NO2",
        "K": "C6H14N2O2",
        "M": "C5H11NO2S",
        "F": "C9H11NO2",
        "T": "C4H9NO3",
        "W": "C11H12N2O2",
        "V": "C5H11NO2",
        "R": "C6H14N4O2",
        "H": "C6H9N3O2",
        "A": "C3H7NO2",
        "N": "C4H8N2O3",
        "D": "C4H7NO4",
        "C": "C3H7NO2S",
        "E": "C5H9NO4",
        "Q": "C5H10N2O3",
        "G": "C2H5NO2",
        "P": "C5H9NO2",
        "S": "C3H7NO3",
        "Y": "C9H11NO3"
    }
    # masses taken from https://www2.chemistry.msu.edu/faculty/reusch/OrgPage/mass.htm
    mass_monoisotopic = {'C': 12.0000, 'H': 1.00783, 'O': 15.9949, 'N': 14.0031, 'S': 31.9721}

    aa_weight_dict = {}
    for aa, formula in amino_acids.items():
        aa_weight = 0
        for atom, num in refindall(r'([A-Z][a-z]*)(\d*)', formula):
            if num == '':
                num = 1
            aa_weight += int(num) * mass_monoisotopic.get(atom)

        aa_weight_dict[aa] = aa_weight
    return aa_weight_dict

# calculate monoisotopic and average mass of peptide
def calculate_masses(peptide):
    aa_weight_dict = prepare_masses_dict()
    # calculate average mass with 'averagine' formula from lecture slides
    # C:4.94 H:7.76 N:1.36 O:1.48 S:0.04
    averagine_mass = (12.0000 * 4.94) + (1.00783 * 7.76) + (14.0031 * 1.36) + (15.9949 * 1.48) + (31.9721 * 0.04)
    # delete mass of H20
    H2O_mass = (len(peptide) - 1) * (1.00783 * 2 + 15.9949)
    monoisotopic_mass = sum(aa_weight_dict.get(aa) for aa in peptide) - H2O_mass
    average_mass = len(peptide) * averagine_mass - H2O_mass
    return monoisotopic_mass, average_mass


def plot_intense_peaks(monoiso_mass):
    # from the slides
    p = [0.55, 0.30, 0.10]
    k = [0, 1, 2]
    x = np.array([monoiso_mass + i for i in k])

    plt.figure(figsize=(10,5))
    plt.bar(x, p, width=0.03)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title("3 most intense peaks")
    plt.savefig("intense_peaks.png")
    plt.show()



def calculate_fragment_ions(peptide):
    b_ions_masses = []
    y_ions_masses = []

    for i in range(1, len(peptide)):
        # delete one O and one H
        OH_mass = 1.00783 + 15.9949
        b_ion_mass = calculate_masses(peptide[:i])[0] - OH_mass
        b_ions_masses.append(b_ion_mass)

        # add one H
        H_mass = 1.00783
        y_ion_mass = calculate_masses(peptide[i:])[0] + H_mass
        y_ions_masses.append(y_ion_mass)

    return b_ions_masses, y_ions_masses


def main(file):
    # task 2
    # read fasta file
    with open(file) as fasta_file:
        sequence = SeqIO.read(fasta_file, "fasta").seq

    peptide_list = trypsin_digestion(sequence)

    peptide_weight_dict = {}
    for peptide in peptide_list:
        peptide_weight_dict[peptide] = calculate_masses(peptide)
    df_peptide_mass = pd.DataFrame.from_dict(peptide_weight_dict, orient='index',
                                             columns=['Monoisotopic Mass', 'Average Mass'])
    df_peptide_mass.index.name = 'Peptide'
    sorted_df_peptide_mass = df_peptide_mass.sort_index(key=lambda x: x.str.len(), ascending=False)
    df_peptide_mass.to_csv('peptide_mass.csv')
    sorted_df_peptide_mass.to_csv('sorted_peptide_mass.csv')
    print('Task 2\n',df_peptide_mass)
    print()
    print(sorted_df_peptide_mass)

    # task 3.1
    peptide_seq = 'HFEEDMGRK'
    monoiso_mass = calculate_masses(peptide_seq)[0]
    print(f'Task 3\nMonoisotopic Mass of {peptide_seq}: {monoiso_mass}')

    # task 3.2
    plot_intense_peaks(monoiso_mass)

    # task 3.3
    b_ions_mass, y_ions_mass = calculate_fragment_ions(peptide_seq)
    b_ions_mass = [f'{ion:.2f}' for ion in b_ions_mass]
    y_ions_mass = [f'{ion:.2f}' for ion in y_ions_mass]
    print('b ions mass: ', '\t'.join(b_ions_mass), '\ny ions mass: ', '\t'.join(y_ions_mass))


if __name__ == '__main__':
    args = create_parser()
    main(args.input)
