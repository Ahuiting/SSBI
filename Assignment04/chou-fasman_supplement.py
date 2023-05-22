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

helix_result = []
i=0
while i < len(seq) - 6:
    w_helix = sum([w_a.get(aa) for aa in seq[i:i + 6]])
    if w_helix >= 4:
        start = i
        print(start)
        end = i + 6

        # left
        for j in range(1, i):
            p_helix = sum([p_a.get(aa) for aa in seq[i - j:i - j + 4]])
            if p_helix < 4:
                start = i - j + 1
                break
        # right
        for k in range(i + 3, len(seq) - 4):
            p_helix = sum([p_a.get(aa) for aa in seq[k:k + 4]])
            if p_helix < 4:
                end = k + 3
                break
        if (start, end) not in helix_result:
            print(start)
            helix_result.append((start, end))
        i = end
    i+=1
print(helix_result)

sheet_result = []
i=0
while i <len(seq) - 5:
    start = i
    end = i + 5
    hb_count = 0
    bb_count = 0
    for aa in seq[i:i + 5]:
        hb_count += 1 if w_b.get(aa) == 2 else 0
        bb_count += 1 if w_b.get(aa) == -1 else 0
    if hb_count >= 3 and bb_count <= 1:
        # left
        for j in range(1, i):
            p_sheet = sum([p_b.get(aa) for aa in seq[i - j:i - j + 4]])
            if p_sheet < 4:
                start = i - j + 1
                break
        # right
        for k in range(i + 3, len(seq) - 4):
            p_sheet = sum([p_b.get(aa) for aa in seq[k:k + 4]])
            if p_sheet < 4:
                end = k + 3
                break
        if (start, end) not in sheet_result:
            sheet_result.append((start, end))
        i=end
    i+=1

print(sheet_result)

sheet_result_list = set(num for start, end in sheet_result for num in range(start, end + 1))
helix_result_list = set(num for start, end in helix_result for num in range(start, end + 1))
overlap = list(helix_result_list.intersection(sheet_result_list))
print(overlap)


def conflict(overlap):
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
        p_a_avg = sum([p_a.get(aa) for aa in seq[end - start + 1]]) / (end - start + 1)
        p_b_avg = sum([p_b.get(aa) for aa in seq[end - start + 1]]) / (end - start + 1)
        conflict_dist[(start, end)] = 'H' if p_a_avg > p_b_avg else 'S'

    return conflict_dist


pred_seq = ['-' for i in range(len(seq))]
for i in range(len(seq)):
    if i in helix_result_list and i in sheet_result_list:
        for k, v in conflict(overlap).items():
            if i in range(k[0], k[1] + 1):
                pred_seq[i] = v
    elif i in helix_result_list:
        pred_seq[i] = 'H'

    elif i in sheet_result_list:
        pred_seq[i] = 'S'

print(sheet_result_list)
print(helix_result_list)
print(''.join(pred_seq))


print('precent',sum(1 for i in range(len(seq)) if pred_seq[i] == ss_ref[i]) * 100 / len(seq))