from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
import argparse

def matrixbuilder(file,nr=0):
    acid_names=[]
    p = PDBParser()
    structure = p.get_structure("1MOT", file)
    model = structure[nr]
    list_of_acids = []
    for chain in model:
        for acid in chain:
            if "N" in acid and "C" in acid and "H" in acid and "O" in acid:
                acid_names.append(acid.get_resname())
                list_of_acids.append(
                    {"N": acid["N"].coord, "C": acid["C"].coord, "H": acid["H"].coord, "O": acid["O"].coord, })

    mtx = np.empty((len(list_of_acids), len(list_of_acids)))
    mtt = np.zeros((len(list_of_acids), len(list_of_acids))) == 1

    for a1 in range(len(list_of_acids)):
        for a2 in range(len(list_of_acids)):
            mtx[a1, a2] = valbetween(list_of_acids[a1], list_of_acids[a2])
            if mtx[a1, a2] < -2.4:
                mtt[a1, a2] = True
    return mtx,mtt,acid_names


def distance(a, b):
    return sum((a - b) ** 2) ** (1 / 2)


def makecsv(names, matrix, filename):
    with open(filename, "w") as f:
        f.write("\t".join(["-"] + names) + "\n")
        for i in range(matrix.shape[0]):
            f.write(names[i] + "\t" + "\t".join([str(x) for x in matrix[i, :]]) + "\n")


def valbetween(acid1, acid2):
    on = 1 / distance(acid1["O"], acid2["N"])
    ch = 1 / distance(acid1["C"], acid2["H"])
    oh = 1 / distance(acid1["O"], acid2["H"])
    cn = 1 / distance(acid1["C"], acid2["N"])
    faktor = 27.888 * 4.1868
    return (on + ch - oh - cn) * faktor


def check(m, a, b):
    try:
        return m[a, b]
    except:
        return False


def checkbridge(m, i, p):
    for j in range(m.shape[0]):
        if i - p == j and j == i + p:
            continue
        if check(m, i - p, j) and check(m, j, i + p):
            return True
    return False


def check_sheet(matr, r, cs):
    for k in range(len(cs)):
        if cs[k] != -1 or not checkbridge(matr, k, r):
            continue
        if checkbridge(matr, k + 1, r) and cs[k + 1] == -1:
            cs[k] = r
            cs[k + 1] = r
    return cs


def check_helix(matr, r, cs):
    for k in range(len(cs)):
        if cs[k] != -1 or not check(matr, k, k + r):
            continue
        if check(matr, k + 1, r + k + 1) and cs[k + 1] == -1:
            cs[k] = r
            cs[k + 1] = r
    return cs

def make_heatmap(matr,title,file):
    fig, ax = plt.subplots()
    ax.set_title(title,fontsize=10)
    plt.pcolor(matr)
    plt.colorbar()
    plt.savefig(file)
    plt.close()

def generate_structure(mtt):
    r_n = {3: "G", 4: "H", 5: "I", -1: "-", 0: "A", 1: "P"}
    cs = [-1 for _ in range(mtt.shape[0])]
    cs = check_helix(mtt, 4, cs)
    cs = check_helix(mtt, 3, cs)
    cs = check_helix(mtt, 5, cs)
    cs = check_sheet(mtt, 0, cs)
    cs = check_sheet(mtt, 1, cs)
    return "".join(r_n[x] for x in cs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='?', const=1, type=str, default="./5jxv.pdb")
    args = parser.parse_args()
    mtx,mtt,acid_names=matrixbuilder(args.i)
    name=generate_structure(mtt)
    makecsv(acid_names, mtx, "./dssp_matrix.tsv")
    make_heatmap(mtx,name,"./heatmap.png")
    print(name)

if __name__ == "__main__":
    main()

