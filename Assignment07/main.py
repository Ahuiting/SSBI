
#from Bio import SeqIO

import pandas as pd
import re

#fasta_sequence = SeqIO.parse(open("PycharmProjects/ssbi_assignment07/P07327.fasta"),'fasta')
#for fs in fasta_sequence:
 #   name, sequence = fs.id, str(fs.seq)

#tryptic digestion at R and K residues

fsseq='MSTAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKMVAVGICGTDDHVVSGTMVTPLPVILGHEAAGIVESVGEGVTTVKPGDKVIPLAIPQCGKCRICKNPESNYCLKNDVSNPQGTLQDGTSRFTCRRKPIHHFLGISTFSQYTVVDENAVAKIDAASPLEKVCLIGCGFSTGYGSAVNVAKVTPGSTCAVFGLGGVGLSAIMGCKAAGAARIIAVDINKDKFAKAKELGATECINPQDYKKPIQEVLKEMTDGGVDFSFEVIGRLDTMMASLLCCHEACGTSVIVGVPPDSQNLSMNPMLLLTGRTWKGAILGGFKSKECVPKLVADFMAKKFSLDALITHVLPFEKINEGFDLLHSGKSIRTILMF'


print(fsseq)
R_cleavage = fsseq.split('R')
K_cleavage = fsseq.split('K')
print('R_cleavage:', R_cleavage)
print('K_cleavage:', K_cleavage)
RK_cleavage = re.split(r'[KR]', fsseq)
print('RK_cleavage:', RK_cleavage)



resulting_peptides=[]

def tryptic_digest(seq):
    for aa in range(0, len(seq)):
        if seq[aa] == 'R' and seq[aa+1] != 'P':
            resulting_peptides.append(seq[0: aa+1])
           # print('R', resulting_peptides)
        if seq[aa] == 'K'and seq[aa+1] != 'P':
            resulting_peptides.append(seq[0:aa+1])
           # print('K', resulting_peptides)


tryptic_digest(fsseq)
print(resulting_peptides)

df = pd.DataFrame(resulting_peptides)
print(df)

# calculate masses of resulting peptides

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

#masses taken from https://www2.chemistry.msu.edu/faculty/reusch/OrgPage/mass.htm

mass_monoisotopic = {'C': 12.0000, 'H': 1.00783, 'O': 15.9949, 'N': 14.0031, 'S': 31.9721}

#def monoisotopic_peptide_mass(resulting_peptide_list):
 #   for peptide in resulting_peptide_list:
  #      for atom in peptide:




#from internet
from re import findall as refindall

def molecular_weight(molecule):
    return sum(mass_monoisotopic[atom] * int(num or '1')
               for atom, num in refindall(r'([A-Z][a-z]*)(\d*)', molecule)
               )

def proteine_mass(proteine):
    return sum(molecular_weight(amino_acids[char]) for char in proteine)

print(proteine_mass('M'))
print('peptide_task3:' ,proteine_mass('HFEEDMGRK'))


#calculate average mass with 'averagine' formula from lecture slides
#C4.94H7.76N1.36O1.48S0.04


averagine_mass = (12.0000 * 4.94) + (1.00783 * 7.76) + (14.0031 * 1.36) + (15.9949 * 1.48) + (31.9721 * 0.04)
print(averagine_mass)

#idea:
#take each peptide and calculate length
#take this length times averagine mass
#average mass in Dalton


def avergae_mass(peptide):
    for peptide in resulting_peptides:
        av_mass =len(peptide) * averagine_mass
        print(av_mass)


avergae_mass(resulting_peptides)


# console output for resulting peptides + monoisotopic mass + average mass

# make table of resulting peptide + monoisotopic mass + average mass



