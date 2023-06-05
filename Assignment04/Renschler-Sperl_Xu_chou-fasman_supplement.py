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

# Author: Huiting Xu and Desirée Renschler-Sperl

# return start, end index
def extend_left_right(start, end, p_dist):
    # left extend
    left_move = 1
    while left_move <= start:
        p_helix = sum(p_dist.get(aa) for aa in seq[start - left_move:start - left_move + 4])
        if p_helix < 4:
            start = start - left_move + 1
            break
        else:
            left_move += 1

    # right extend
    right_move = 1
    while end + right_move <= len(seq):
        p_helix = sum(p_dist.get(aa) for aa in seq[end + right_move - 4:end + right_move])
        if p_helix < 4:
            end = end + right_move - 1
            break
        else:
            right_move += 1
    return start, end


# return a dict, key is the range of overlap, value is 'H' or 'S'
def conflict_resolve(overlap):
    overlap_tuple = []
    start_num = overlap[0]
    for i in range(1, len(overlap)):
        if overlap[i] != overlap[i - 1] + 1:
            overlap_tuple.append((start_num, overlap[i - 1]))
            start_num = overlap[i]
        if i == len(overlap) - 1:
            overlap_tuple.append((start_num, overlap[i]))
    conflict_dist = {}
    for start, end in overlap_tuple:
        # average calculate
        p_a_avg = sum(p_a.get(aa) for aa in seq[end - start + 1]) / (end - start + 1)
        p_b_avg = sum(p_b.get(aa) for aa in seq[end - start + 1]) / (end - start + 1)
        conflict_dist[(start, end)] = 'H' if p_a_avg > p_b_avg else 'S'
    return conflict_dist

if __name__ == '__main__':
    # helix
    helix_cores = []
    idx_a = 0
    # find the helix core
    while idx_a < len(seq) - 6:
        w_helix = sum(w_a.get(aa) for aa in seq[idx_a:idx_a + 6])
        if w_helix >= 4:
            start, end = extend_left_right(idx_a, idx_a + 6, p_a)
            helix_cores.extend(list(range(start, end + 1)))
            idx_a = idx_a + 6 # skip already found and extended ones while searching for cores
        idx_a += 1
    print('Helix cores: ',''.join(seq[i] if i in helix_cores else '-'for i in range(len(seq))))
    print('Helix: ',''.join('H' if i in helix_cores else '-' for i in range(len(seq))))

    # sheet
    sheet_cores = []
    idx_b = 0
    # find the sheet core
    while idx_b < len(seq) - 5:
        start = idx_b
        end = idx_b + 5

        hb_count = [w_b.get(aa) for aa in seq[start:end]].count(2)
        bb_count = [w_b.get(aa) for aa in seq[start:end]].count(-1)

        if hb_count >= 3 and bb_count <= 1:
            start, end = extend_left_right(start, end, p_b)
            sheet_cores.extend(list(range(start, end + 1)))
            idx_b = idx_b + 5 # skip already found and extended ones while searching for cores
        idx_b += 1
    print('Sheet cores: ', ''.join(seq[i] if i in sheet_cores else '-' for i in range(len(seq))))
    print('Sheet: ',''.join('S' if i in sheet_cores else '-' for i in range(len(seq))))

    overlap_list = [element for element in helix_cores if element in sheet_cores]

    pred_seq = ['-' for i in range(len(seq))]
    for idx_a in range(len(seq)):
        if idx_a in overlap_list:
            for k, v in conflict_resolve(overlap_list).items():
                if idx_a in range(k[0], k[1] + 1):
                    pred_seq[idx_a] = v
        elif idx_a in helix_cores:
            pred_seq[idx_a] = 'H'
        elif idx_a in sheet_cores:
            pred_seq[idx_a] = 'S'


    print('Prediction: ',''.join(pred_seq))
    print('Reference: ', ''.join(pred_seq))
    print('Accuracy: ',sum(1 if pred_seq[i] == ss_ref[i] else 0 for i in range(len(seq))) / len(seq)* 100,'%')

