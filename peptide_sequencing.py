#from itertools import product


def get_codon_table():
    d = {}
    d['CAU'] = 'H'; d['CAC'] = 'H'; d['CAA'] = 'Q'; d['CAG'] = 'Q'; d['CCU'] = 'P'; d['CCC'] = 'P'; d['CCA'] = 'P'; d['CCG'] = 'P';
    d['CGU'] = 'R'; d['CGC'] = 'R'; d['CGA'] = 'R'; d['CGG'] = 'R'; d['CUU'] = 'L'; d['CUC'] = 'L'; d['CUA'] = 'L'; d['CUG'] = 'L';
    d['GAU'] = 'D'; d['GAC'] = 'D'; d['GAA'] = 'E'; d['GAG'] = 'E'; d['GCU'] = 'A'; d['GCC'] = 'A'; d['GCA'] = 'A'; d['GCG'] = 'A';
    d['GGU'] = 'G'; d['GGC'] = 'G'; d['GGA'] = 'G'; d['GGG'] = 'G'; d['GUU'] = 'V'; d['GUC'] = 'V'; d['GUA'] = 'V'; d['GUG'] = 'V';
    d['UAU'] = 'Y'; d['UAC'] = 'Y'; d['UAA'] = '*'; d['UAG'] = '*'; d['UCU'] = 'S'; d['UCC'] = 'S'; d['UCA'] = 'S'; d['UCG'] = 'S';
    d['UGU'] = 'C'; d['UGC'] = 'C'; d['UGA'] = '*'; d['UGG'] = 'W'; d['UUU'] = 'F'; d['UUC'] = 'F'; d['UUA'] = 'L'; d['UUG'] = 'L';
    d['AAU'] = 'N'; d['AAC'] = 'N'; d['AAA'] = 'K'; d['AAG'] = 'K'; d['ACU'] = 'T'; d['ACC'] = 'T'; d['ACA'] = 'T'; d['ACG'] = 'T';
    d['AGU'] = 'S'; d['AGC'] = 'S'; d['AGA'] = 'R'; d['AGG'] = 'R'; d['AUU'] = 'I'; d['AUC'] = 'I'; d['AUA'] = 'I'; d['AUG'] = 'M';
    return d


'''def get_amino_table():
    d = {}
    d['H'] = ['CAU', 'CAC']; d['Q'] = ['CAA', 'CAG']; d['P'] = ['CCU', 'CCC', 'CCA', 'CCG']; d['D'] = ['GAU', 'GAC']
    d['R'] = ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']; d['L'] = ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG']
    d['E'] = ['GAA', 'GAG']; d['A'] = ['GCU', 'GCC', 'GCA', 'GCG']; d['G'] = ['GGU', 'GGC', 'GGA', 'GGG']
    d['V'] = ['GUU', 'GUC', 'GUA', 'GUG']; d['Y'] = ('UAU', 'UAC'); d['*'] = ['UAA', 'UAG', 'UGA']
    d['S'] = ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']; d['C'] = ['UGU', 'UGC']; d['W'] = ['UGG']
    d['F'] = ['UUU', 'UUC']; d['N'] = ['AAU', 'AAC']; d['K'] = ['AAA', 'AAG']; d['T'] = ['ACU', 'ACC', 'ACA', 'ACG']
    d['I'] = ['AUU', 'AUC']; d['M'] = ['AUG']
    return d'''


def protein_translation(rna_string):
    amino_acid = ''
    codon_table = get_codon_table()
    i = 0
    while i < len(rna_string):
        current = rna_string[i: i+3]
        if codon_table[current] == '*':
            break
        amino_acid += codon_table[current]
        i += 3
    return amino_acid


def dna_to_rna(dna):
    rna = ''
    for nucleotide in dna:
        if nucleotide == 'T':
            rna += 'U'
        else:
            rna += nucleotide
    return rna


def reverse_rna(rna):
    res = ''
    for i in range(len(rna) - 1, -1, -1):
        if rna[i] == 'A':
            res += 'U'
        elif rna[i] == 'U':
            res += 'A'
        elif rna[i] == 'G':
            res += 'C'
        elif rna[i] == 'C':
            res += 'G'
    return res


'''def peptide_encoding_problem1(dna, amino_acid):
    amino_table = get_amino_table()
    codons_in_amino_acid = []
    for el in amino_acid:
        codons_in_amino_acid.append(amino_table[el])

    transcriptions = list(product(*codons_in_amino_acid))
    for i in range(len(transcriptions)):
        transcriptions[i] = ''.join(transcriptions[i])

    result = []
    rna = dna_to_rna(dna)
    size = len(amino_acid) * 3
    length = len(dna)
    for i in range(length-size+1):
        current = rna[i: i+size]
        if current in transcriptions:
            result.append(dna[i: i+size])

    numbers = []
    complement = reverse_rna(rna)
    for i in range(length-size+1):
        current = complement[i: i+size]
        if current in transcriptions:
            result.append(dna[length-i-size: length-i])

    return result'''


def peptide_encoding(dna, amino_acid):
    rna = dna_to_rna(dna)
    window_size = len(amino_acid)*3
    result = []
    for i in range(len(dna)-window_size+1):
        if protein_translation(rna[i: i+window_size]) == amino_acid:
            result.append(dna[i: i+window_size])
        if protein_translation(reverse_rna(rna[i: i+window_size])) == amino_acid:
            result.append(dna[i: i + window_size])
    return result


def subpeptides_count(length):
    return length * (length - 1)


def get_mass_table():
    return {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
            'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


def get_mass(amino_acid):
    summ = 0
    mass_table = get_mass_table()
    for el in amino_acid:
        summ += mass_table[el]
    return summ


def get_theoretical_spectrum(amino_acid):
    spectrum = [0]
    length = len(amino_acid)
    for window_size in range(1, length):
        for i in range(1, length+1):
            if i-window_size < 0:
                spectrum.append(get_mass(amino_acid[i-window_size:] + amino_acid[: i]))
            else:
                spectrum.append(get_mass(amino_acid[i - window_size: i]))
    spectrum.append(get_mass(amino_acid))
    spectrum.sort()
    return spectrum


def peptide_count(mass):
    masses = (57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186)
    N = {0: 1}
    for i in range(57, 243):
        N[i] = 0
        for j in range(len(masses)):
            k = i - masses[j]
            if k in N:
                N[i] += N[k]
    for i in range(243, mass+1):
        N[i] = 0
        for j in range(len(masses)):
            N[i] += N[i-masses[j]]
    return N[mass]


def main():
    rna_string = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    print(protein_translation(rna_string))

    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    amino_acid = 'MA'
    print(peptide_encoding(dna, amino_acid))

    length = 34215
    print(subpeptides_count(length))

    amino_acid = 'LEQN'
    print(get_theoretical_spectrum(amino_acid))

    mass = 1024
    print(peptide_count(mass))


main()
