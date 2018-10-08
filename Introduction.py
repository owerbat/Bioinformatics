NUCLEOTIDES = ('A', 'T', 'G', 'C')


def pattern_count(pattern, genome):
    res = 0
    for i in range(len(genome)-len(pattern)):
        if pattern == genome[i: i+len(pattern)]:
            res += 1
    return res


def frequent_words(k, genome):
    frequent_list = []
    dic = {}
    for i in range(len(genome)-k):
        current = genome[i: i+k]
        if current not in dic:
            dic.update({current: pattern_count(current, genome)})
    maximum = max([dic[key] for key in dic])
    for key in dic:
        if dic[key] == maximum:
            frequent_list.append(key)
    return frequent_list


def reverse_complement(genome):
    res = ''
    for i in range(len(genome)-1, -1, -1):
        if genome[i] == 'A':
            res += 'T'
        elif genome[i] == 'T':
            res += 'A'
        elif genome[i] == 'G':
            res += 'C'
        elif genome[i] == 'C':
            res += 'G'
    return res


def main():
    pattern = 'ATAT'
    genome = 'GATATATGCATATACTT'
    print(pattern_count(pattern, genome))

    k = 4
    genome2 = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    print(frequent_words(k, genome2))

    genome3 = 'AAAACCCGGT'
    print(reverse_complement(genome3))


main()
