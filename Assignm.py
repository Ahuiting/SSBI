# Author: Huiting Xu
import argparse

import pandas as pd
from Bio import SeqIO
import numpy as np


def create_parser():
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', help='PATH_TO_FASTA_FILE')
    p.add_argument('--min-loop-length', default=3, help='MIN-LOOP-LENGTH')
    p.add_argument('--score-GC', default=3, help='SCORE-GC')
    p.add_argument('--score-AU', default=2, help='SCORE-AU')
    p.add_argument('--score-GU', default=1, help='SCORE-GU')


def delta_score(tup):
    if tup in [('G', 'C'), ('C', 'G')]:
        return 1
    elif tup in [('A', 'U'), ('U', 'A')]:
        return 1
    elif tup in [('G', 'U'), ('U', 'G')]:
        return 1
    else:
        return 0


def DPmatrix(seq, min_loop_size):
    # initialize
    matrix = np.zeros((len(seq), len(seq)))
    # 4 cases
    for l in range(min_loop_size + 1, len(seq)):
        for j in range(l, len(seq)):
            i = j - l;
            case1 = matrix[i + 1, j];
            case2 = matrix[i, j - 1];
            case3 = matrix[i + 1, j - 1] + delta_score((seq[i], seq[j]));
            if j - i > min_loop_size + 1:
                case4_scores = []
                for k in range(i + 1 + min_loop_size, j):
                    case4_scores.append(matrix[i, k] + matrix[k + 1, j]);
                case4 = max(case4_scores);
                matrix[i, j] = max(case1, case2, case3, case4);
            else:
                matrix[i, j] = max(case1, case2, case3);

    return matrix;





if __name__ == '__main__':
    args = create_parser()
    seq=''
    for seq_record in SeqIO.parse('test.fasta', "fasta"):
        seq = str(seq_record.seq)
    matrix=DPmatrix(seq,2)
    df = pd.DataFrame(matrix, index=list(seq), columns=list(seq))
    print(df)