# ATAT
# GATATATGCATATACTT


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
    for i in NUCLEOTIDES:
        for j in NUCLEOTIDES:
            for k in NUCLEOTIDES:
                for l in NUCLEOTIDES:
                    string = i+j+k+l
                    dic.update({string: pattern_count(string, genome)})
    maximum = max([dic[key] for key in dic])
    for key in dic:
        if dic[key] == maximum:
            frequent_list.append(key)
    return frequent_list


def main():
    pattern = 'ATAT'
    genome = 'GATATATGCATATACTT'
    print(pattern_count(pattern, genome))

    k = 4
    genome2 = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    lst = frequent_words(k, genome2)
    lst.sort()
    print(lst[0], lst[1])



main()
