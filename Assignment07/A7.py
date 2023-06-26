# Author: Huiting Xu and Desirée Renschler-Sperl
from Bio import SeqIO


def trypsin_digestion(sequence):
    peptides = []
    peptide = ""
    previous_aa = ""

    for aa in sequence:
        if (aa == "K" or aa == "R") and previous_aa != "P":
            peptides.append(peptide + aa)
            peptide = ""
        else:
            peptide += aa
        previous_aa = aa
    if peptide:
        peptides.append(peptide)
    return peptides


def main(file):
    with open(file) as fasta_file:
        sequence = SeqIO.read(fasta_file, "fasta").seq
    print(trypsin_digestion(sequence))


if __name__ == '__main__':
    file = "P07327.fasta"
    main(file)
