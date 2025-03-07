# Author: Huiting Xu and Desirée Renschler-Sper
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
                    traceback(i, k)
                    traceback(k + 1, j)


def print_bracket_dot(seq, pairs):
    output = ["." for i in range(len(seq))]
    for pair in pairs:
        list(pair).sort()
        output[pair[0]] = '('
        output[pair[1]] = ')'
    return ''.join(output)


def print_bqseq(seq, pairs):
    output = ['\t'.join([str(i + 1), seq[i], str(0)]) for i in range(len(seq))]
    for pair in pairs:
        list(pair).sort()
        output[pair[0]] = '\t'.join([str(pair[0] + 1), seq[pair[0]], str(pair[1] + 1)])
        output[pair[1]] = '\t'.join([str(pair[1] + 1), seq[pair[1]], str(pair[0] + 1)])
    return '\n'.join(output)


if __name__ == '__main__':
    'try -i test.fasta --min-loop-length 5 --score-GC 10 --score-AU 5 --score-GU 2'

    args = create_parser()
    score_GC = args.score_GC
    score_AU = args.score_AU
    score_GU = args.score_GU
    min_loop_len = args.min_loop_length
    file = args.input

    seq = ''
    for seq_record in SeqIO.parse(file, "fasta"):
        seq = str(seq_record.seq)
    matrix = DPmatrix(seq)
    df = pd.DataFrame(matrix, index=list(seq), columns=list(seq))
    print(df)

    pairs = []
    traceback(0, len(seq) - 1)
    print(pairs)

    print(print_bracket_dot(seq, pairs))
    print(print_bqseq(seq, pairs))

    with open(file[:-6] + ".bpseq", "w") as f:
        f.write("Filename: " + file + "\nMin-loop: " + str(min_loop_len) + "\nGC: " + str(score_GC) + "\nAU: " + str(
            score_AU) + "\nGU: " + str(score_GU) + "\nScore: " + str(matrix[0][len(seq) - 1]))
        f.write("\n" + print_bqseq(seq, pairs))
        f.write("\n" + print_bracket_dot(seq, pairs))
