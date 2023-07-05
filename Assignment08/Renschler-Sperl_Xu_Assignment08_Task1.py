# Author: Huiting Xu and Desirée Renschler-Sperl

import argparse
import numpy as np


# create perse
def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', help='txt_file')
    return p.parse_args()


# read file and return a list of numbers
def read_file(file_path):
    numbers = []
    with open(file_path, 'r') as file:
        for line in file:
            number = float(line.strip())
            numbers.append(number)
    return numbers


# check if the number is in the range of mass of amino acid,  if true, return the value, else False
def check_range(list, number):
    aa = []
    for i in range(0, len(list)):
        if abs(list[i] - number) <= 0.055 * 2:
            #return list[i]
            aa.append(list[i])
    if len(aa) != 0:
        return aa
    else:
        return False


# return the b_ions, y_ions indices of ion_list
def separate_spectrum(ion_list, mono_iso_dict):
    idx_list = []
    # get all possible cases that the difference of two values is mass of a amino acid
    for i in range(len(ion_list)):
        for j in range(i + 1, len(ion_list)):
            if check_range(list(mono_iso_dict.keys()), ion_list[j] - ion_list[i]):
                idx_list.append([i, j])
    # separate numbers to two part
    idx_arr = np.array(idx_list)
    ion_1 = idx_list[0]
    ion_2 = idx_list[1]
    for i in range(2, len(idx_arr)):
        a, b = idx_arr[i]
        if np.count_nonzero(idx_arr[:, 0] == a) == 1 \
                or np.count_nonzero(idx_arr[:, 1] == b) == 1:
            if ion_1[-1] == a:
                ion_1.append(b)
            if ion_2[-1] == a:
                ion_2.append(b)
    # check which list is b ions which is y ions
    if check_range(list(mono_iso_dict.keys()), ion_list[ion_1[0]] - 1):
        y_ions = ion_1
        b_ions = ion_2
    else:
        y_ions = ion_2
        b_ions = ion_1
    return b_ions, y_ions


# from index map to amino acid
def calculate_peptide(sublist, ion_list, mono_iso_dict):
    peptides = []
    for i in range(len(sublist) - 1):
        k = check_range(list(mono_iso_dict.keys()), ion_list[sublist[i + 1]] - ion_list[sublist[i]])
        if len(k)==1:
            peptides.append(mono_iso_dict.get(k[0]))
        else:
            aa=[]
            for i in k:
                aa.extend(mono_iso_dict.get(i))
            peptides.append(aa)

    return peptides


# get all the possible peptide
def output_peptide(peptide_list):
    result = ['']
    for sublist in peptide_list:
        result = [x for x in result for _ in range(len(sublist))]
        for i in range(len(result)):
            if len(sublist) == 1:
                result[i] += sublist[0]
            else:
                result[i] += sublist[i % 2]  # for alternative cases
    return result


def main(filepath):
    mono_iso_dict = {71.03711: ['A'], 156.10111: ['R'], 114.04293: ['N'], 115.02694: ['D'],
                     103.00919: ['C'], 129.04259: ['E'], 128.05858: ['Q'], 57.02146: ['G'],
                     137.05891: ['H'], 113.08406: ['I', 'L'], 128.09496: ['K'], 131.04049: ['M'],
                     147.06841: ['F'], 97.05276: ['P'], 87.03203: ['S'], 101.04768: ['T'],
                     186.07931: ['W'], 163.06333: ['Y'], 99.06841: ['V']}

    ion_list = (read_file(filepath))
    b_ion_list, y_ion_list = separate_spectrum(ion_list, mono_iso_dict)
    aa_list = calculate_peptide(b_ion_list, ion_list, mono_iso_dict)
    # add the last aa in y_ions as first aa of peptide
    aa_list.insert(0, calculate_peptide(y_ion_list, ion_list, mono_iso_dict)[-1:][0])
    all_possible_cases = output_peptide(aa_list)
    print('There are', len(all_possible_cases), 'possible peptides:')
    print('\n'.join(all_possible_cases))
    with open('all_possible_peptides.txt', 'w') as f:
        f.write('\n'.join(all_possible_cases))



if __name__ == '__main__':
    args = create_parser()
    main(args.input)
