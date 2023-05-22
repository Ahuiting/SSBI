# helix p_alpha relative probabilities
p_a = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12,  # builders
    'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77,  # indifferent
    'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53  # breakers
}

# helix builder/indifferent/breaker class weights
w_a = {aa: (1 if score > 1.1 else -1 if score < 0.75 else 0.5) for aa, score in p_a.items()}
# strand p_beta relative probabilities
p_b = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19,
    # builders
    'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80,  # indifferent
    'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26  # breakers
}

# strand builder/indifferent/breaker class weights
# changed to 2/0/-1 encoding for convenience, so that a window of 5 is a strand core if sum(window) >= 5.
w_b = {aa: (2 if score > 1 else -1 if score < 0.78 else 0) for aa, score in p_b.items()}

# sequence of 5JXV
seq = 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'

# reference secondary structure for accuracy calculation
ss_ref = '-SSSSSSS----SSSSSSS---HHHHHHHHHHHHHH-----SSSSS----SSSSS-'


def getarea(seq, core_tbl, extend_p, window_core, window_extend, thresh_core, thresh_exp, letter):
    regions = ["-" for _ in seq]
    idx = 0
    while idx < len(seq) - window_core:
        core_score = sum(core_tbl[i] for i in seq[idx:idx + window_core])
        if core_score >= thresh_core:
            left = 1
            right = 1
            while idx - left > 0 and thresh_exp < sum(
                    extend_p[i] for i in seq[idx - left:idx - left + window_extend]) / window_extend:
                left += 1
            while idx + right + window_core < len(seq) and \
                    thresh_exp < sum(extend_p[i] for i in seq[idx + right + window_core:window_core + idx + right + window_extend]) / window_extend:
                right += 1
            for i in range(idx - left + 1, idx + right + window_core):
                regions[i] = letter
            idx += 1
        else:
            idx += 1
    return "".join(regions)


def conflict(r_a, r_b, p_a, p_b, seq):
    combi = []
    idx = 0
    while idx < len(r_a):
        if r_a[idx] != "-" and r_b[idx] != "-":
            toright = 1
            while idx + toright < len(r_a) and r_a[idx + toright] != "-" and r_b[idx + toright] != "-":
                toright += 1
            winner = r_b[idx]
            if sum(p_a[i] for i in seq[idx:idx + toright]) > sum(p_b[i] for i in seq[idx:idx + toright]):
                winner = r_a[idx]
            for i in range(idx, idx + toright):
                combi.append(winner)
            idx += toright
        elif r_a[idx] == "-":
            combi.append(r_b[idx])
            idx += 1
        else:
            combi.append(r_a[idx])
            idx += 1
    return "".join(combi)


helix = getarea(seq, w_a, p_a, 6, 4, 4, 1, "H")
sheet = getarea(seq, w_b, p_b, 5, 4, 4, 1, "S")
cmb = conflict(helix, sheet, p_a, p_b, seq)
q3 = sum(1 for i in range(len(seq)) if cmb[i] == ss_ref[i]) * 100 / len(seq)
print(seq)
print(helix)
print(sheet)
print(cmb)
print(ss_ref)
print("Q3:", q3, "%")
exit(0)
