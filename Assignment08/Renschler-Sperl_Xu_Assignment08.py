# Author: Huiting Xu and Desirée Renschler-Sperl
import argparse

import numpy as np


def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', help='txt_file')
    return p.parse_args()


def read_file(file_path):
    numbers = []
    with open(file_path, 'r') as file:
        for line in file:
            number = float(line.strip())
            numbers.append(number)
    return numbers


def calculate_peptide_sequences(mz_values):
    peptide_sequences = []
    num_ions = len(mz_values)
    for i in range(num_ions):
        for j in range(i + 1, num_ions):
            b_ion_mass = mz_values[i]
            y_ion_mass = mz_values[j]
            peptide_sequence = f'{b_ion_mass}+{y_ion_mass}'
            peptide_sequences.append(peptide_sequence)
    return peptide_sequences


def identify_peptide(spectrum):
    peptide_sequence = ''
    list = []
    for i in range(len(spectrum) - 1):
        ion1 = spectrum[i]
        ion2 = spectrum[i + 1]
        mass_diff = ion2 - ion1
        list.append(mass_diff)
    return list


def check_numbers(list, number):
    for i in range(0, len(list)):
        if abs(list[i] - number) <= 0.055 * 2:
            return True
    return False


if __name__ == '__main__':

    mono_isotopic_masses = np.array(
        [71.03711, 156.10111, 114.04293, 115.02694, 103.00919, 129.04259, 128.05858, 57.02146,
         137.05891, 113.08406, 113.08406, 128.09496, 131.04049, 147.06841, 97.05276, 87.03203,
         101.04768, 186.07931, 163.06333, 99.06841])
    ion_list = (read_file('materials/b_y_spectrum.txt'))

    idx_list = []
    for i in range(len(ion_list)):
        for j in range(i + 1, len(ion_list)):
            if (check_numbers(mono_isotopic_masses, ion_list[j] - ion_list[i])):
                idx_list.append([i, j])
                # all possible cases
    idx_arr = arr = np.array(idx_list)
    ion_1 = idx_list[0]
    ion_2 = idx_list[1]
    for i in range(2, len(idx_arr)):
        a, b = idx_arr[i]
        if np.count_nonzero(idx_arr[:, 0] == a) == 1 or np.count_nonzero(
                idx_arr[:, 1] == b) == 1:
            if ion_1[-1] == a:
                ion_1.append(b)
            if ion_2[-1] == a:
                ion_2.append(b)

    print(ion_2)
    print(ion_1)
