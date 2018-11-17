import peptide_sequencing as ps


def cyclopeptide_score(peptide, experimental_spectrum):
    theoretical_spectrum = ps.get_theoretical_spectrum(peptide)
    score = 0
    for value in experimental_spectrum:
        if value in theoretical_spectrum:
            score += 1
            theoretical_spectrum.remove(value)
    return score


def mass(peptide):
    summ = 0
    mass_table = ps.get_mass_table()
    for el in peptide:
        summ += mass_table[el]
    return summ


def parent_mass(spectrum):
    return max(spectrum)


def expand(peptides):
    masses = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    if peptides != ['*']:
        result = []
        for peptide in peptides:
            for value in masses:
                result.append(peptide+value)
        return result
    else:
        return masses


def get_linear_spectrum(peptide):
    spectrum = [0]
    for window_size in range(1, len(peptide)+1):
        for i in range(len(peptide)-window_size+1):
            spectrum.append(mass(peptide[i: i+window_size]))
    spectrum.sort()
    return spectrum


def get_cyclospectrum(peptide):
    spectrum = [0]
    length = len(peptide)
    for window_size in range(1, length):
        for i in range(1, length + 1):
            if i - window_size < 0:
                spectrum.append(mass(peptide[i - window_size:] + peptide[: i]))
            else:
                spectrum.append(mass(peptide[i - window_size: i]))
    spectrum.append(mass(peptide))
    spectrum.sort()
    return spectrum


def is_consistent(peptide, spectrum):
    theoretical_spectrum = get_linear_spectrum(peptide)
    experimental_spectrum = spectrum[:]
    for value in theoretical_spectrum:
        if value in experimental_spectrum:
            experimental_spectrum.remove(value)
        else:
            return False
    return True


def cyclopeptide_sequencing(spectrum):
    peptides = ['*']
    result = []
    while len(peptides) > 0:
        peptides = expand(peptides)
        to_delete = []
        for peptide in peptides:
            if mass(peptide) == parent_mass(spectrum):
                if get_cyclospectrum(peptide) == spectrum:
                    result.append(peptide)
                to_delete.append(peptide)
            elif not is_consistent(peptide, spectrum):
                to_delete.append(peptide)
        for el in to_delete:
            peptides.remove(el)
    return result


def peptide_to_mass_string(peptide):
    result = ''
    for amino_acid in peptide:
        result += str(mass(amino_acid)) + '-'
    return result[: len(result)-1]


def main():
    peptide = 'NQEL'
    experimental_spectrum = (0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484)
    print(cyclopeptide_score(peptide, experimental_spectrum))

    spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
    peptides = cyclopeptide_sequencing(spectrum)
    lst = list(set([peptide_to_mass_string(el) for el in peptides]))
    print(lst)


main()
