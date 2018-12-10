from itertools import product


def diff_count(pattern1, pattern2):
    if len(pattern1) != len(pattern2):
        return -1

    error_count = 0
    for i, item in enumerate(pattern1):
        if pattern1[i] != pattern2[i]:
            error_count += 1
    return error_count


def get_patterns(k_mer, d):
    patterns = []
    for item in product('ATGC', repeat=len(k_mer)):
        if 0 <= diff_count(item, k_mer) <= d:
            patterns.append(''.join(item))
    return patterns


def motif_enumeration(dna_list, k, d):
    patterns = []
    for dna in dna_list:
        for i in range(len(dna)-k+1):
            k_mer = dna[i: i+k]
            for pattern in get_patterns(k_mer, d):
                count = 0
                for _dna in dna_list:
                    for j in range(len(_dna) - k + 1):
                        if 0 <= diff_count(_dna[j: j + k], pattern) <= d:
                            count += 1
                            break
                if count == len(dna_list):
                    patterns.append(pattern)
    return list(set(patterns))


def main():
    dna_list = ['CACTGATCGACTTATC', 'CTCCGACTTACCCCAC', 'GTCTATCCCTGATGGC', 'CAGGGTTGTCTTGTCT']
    k = 4
    d = 1
    print(motif_enumeration(dna_list, k, d))


main()
