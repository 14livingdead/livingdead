def D(Pattern, Dna):
    k = len(Pattern)
    dist = 0
    for Seq in Dna:
        hd = k
        for j in range(0, len(Seq)-k):
            kmer = Seq[j:j+k]
            newHam = HammingDistanz(Pattern, kmer)
            if hd > newHam:
                hd = newHam
        dist = dist + hd
    return dist

def HammingDistanz(Pattern, kmer):
    hd = 0
    for j, c in enumerate(Pattern):
        if c != kmer[j]:
            hd += 1
    return hd

def NuToS(basenzahl):
    return {
        0: 'A',
        1: 'C',
        2: 'G',
        3: 'T'
    }[basenzahl]

def NTP(basenzahl,k):
    if k == 1:
        return NuToS(basenzahl)
    PrefixBasenzahl = int(basenzahl / 4)
    r = int(basenzahl % 4)
    symbol = NuToS(r)
    PrefixPattern = NTP(PrefixBasenzahl, k-1)
    return PrefixPattern + symbol

def MedString (k, Dna):
    Dna = Dna.split(' ')
    dist = k*len(Dna)
    best_kmer = ""
    for j in range(0, (4**k-1)):
        Pattern = NTP(j, k)
        PattDistance = D(Pattern, Dna)
        if PattDistance <= dist:
            dist = PattDistance
            best_kmer = Pattern
    return best_kmer

print (MedString (6, "TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG"))
