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
    p.add_argument('--min-loop-length', type=int, default=3, help='MIN-LOOP-LENGTH')
    p.add_argument('--score-GC', type=int, default=3, help='SCORE-GC')
    p.add_argument('--score-AU', type=int, default=2, help='SCORE-AU')
    p.add_argument('--score-GU', type=int, default=1, help='SCORE-GU')
    return p.parse_args()


def delta_score(tup):
    if tup in [('G', 'C'), ('C', 'G')]:
        return score_GC
    elif tup in [('A', 'U'), ('U', 'A')]:
        return score_AU
    elif tup in [('G', 'U'), ('U', 'G')]:
        return score_GU
    else:
        return 0


def DPmatrix(seq):
    # initialize
    matrix = np.zeros((len(seq), len(seq)))
    # 4 cases
    for l in range(1, len(seq)):
        for j in range(l, len(seq)):
            i = j - l
            case1 = matrix[i + 1, j]
            case2 = matrix[i, j - 1]
            case3 = matrix[i + 1, j - 1] + delta_score((seq[i], seq[j])) if j - i >= min_loop_len + 1 else 0

            case4_scores = []
            for k in range(i + 1, j):
                case4_scores.append(matrix[i, k] + matrix[k + 1, j])
            case4 = max(case4_scores) if case4_scores else 0

            matrix[i, j] = max(case1, case2, case3, case4)

    return matrix;


def traceback(i, j):
    if i < j:
        if matrix[i][j] == matrix[i + 1, j]:
            traceback(i + 1, j)
        elif matrix[i][j] == matrix[i, j - 1]:
            traceback(i, j - 1)
        elif j - i >= min_loop_len + 1 and matrix[i][j] == matrix[i + 1, j - 1] + delta_score((seq[i], seq[j])):
            pairs.append((i, j))
            traceback(i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if matrix[i][j] == matrix[i, k] + matrix[k + 1, j]:
                    traceback(i.k)
                    traceback(k + 1, j)


if __name__ == '__main__':

    args = create_parser()
    score_GC = args.score_GC
    score_AU = args.score_AU
    score_GU = args.score_GU
    min_loop_len = args.min_loop_length

    seq = ''
    for seq_record in SeqIO.parse(args.input, "fasta"):
        seq = str(seq_record.seq)
    matrix = DPmatrix(seq)
    df = pd.DataFrame(matrix, index=list(seq), columns=list(seq))
    print(df)

    pairs = []
    traceback(0, len(seq) - 1)
    print(pairs)
