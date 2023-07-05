# Author: Huiting Xu and Desirée Renschler-Sperl
import argparse

import pandas as pd

# create perse
def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', help='txt_file')
    return p.parse_args()

def read_csv(filepath):
    with open(filepath, 'r') as file:
        for index, row in enumerate(file):
            if row.startswith('prot_hit_num'):
                break
    df = pd.read_csv('FToArGsTR.csv', skiprows=index, index_col=False)
    return df


def count_species(df):
    species = [i.split('_')[-1] for i in df['prot_acc']]
    species_dict = dict((i, species.count(i)) for i in species)
    return species_dict

def main(path):
    df = read_csv(path)
    print('----Task 2.1----')
    print(pd.DataFrame.from_dict(count_species(df), orient='index', columns=['counts']))
    print('----Task 2.2----')
    print(df['prot_hit_num'].iloc[-1])
    print('----Task 2.3----')
    print(len(df[(df['prot_matches'] >= 2) & (df['prot_family_member'] == 1)]['prot_acc'].unique()))
    print('----Task 2.6----')
    print(len(df[df['prot_family_member'] >= 1]['pep_seq'].unique()))

if __name__ == '__main__':
    args = create_parser()
    main(args.input)


