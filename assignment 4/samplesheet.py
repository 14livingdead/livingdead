import random

def paToNu(pattern):
    if len(pattern) == 0:
        return 0
    return 4 * paToNu(pattern[0:-1]) + SyToNu(pattern[-1:])

def SyToNu(symbol):
    return {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }[symbol]

def NTP(index, k):
    if k == 1:
        return NuToS(index)
    return NTP(index // 4, k-1) + NuToS(index % 4)

def NuToS(index):
    return {
        0: 'A',
        1: 'C',
        2: 'G',
        3: 'T'
    }[index]

def profileProbable(text, k, profile):
    maxprob = 0
    kmer = text[0:k]
    for i in range(0, len(text) - k +1):
        prob =1
        pattern =text[i:i+k]
        for j in range(k):
            l = SyToNu(pattern[j])
            prob *= profile [l][j]
        if prob > maxprob:
            maxprob =prob
            kmer = pattern
    return kmer


def profileForm(motifs):
    k= len(motifs[0])
    profile = [[1 for i in range(k)] for j in range(4)]
    for index in motifs:
        for i in range(len(index)):
            j = SyToNu(index[i])
            profile[j][i] +=1
    for index in profile:
        for i in range(len(index)):
            index[i] = index[i]/len(motifs)
    return profile

def consensus(profile):
    str = ""
    for i in range(len(profile[0])):
        max = 0
        loc = 0
        for j in range(4):
            if profile[j][i] > max:
                loc = j
                max = profile[j][i]
        str+=NuToS(loc)
    return str

def score(motifs):
    profile = profileForm(motifs)
    cons = consensus(profile)
    score = 0
    for index in motifs:
        for i in range(len(index)):
            if cons[i] != index[i]:
                score +=1
    return score

def randomMotifSearch(DNA, k, t):
    bestMotifs = []
    motifs = []
    for index in range(t):
        random.seed()
        i= random.randint(0, len(DNA[index])-k)
        motifs.append(DNA[index][i:i+k])
    bestMotifs = motifs.copy()
    count = 0
    while True:
        profile = profileForm(motifs)
        for index in range(t):
            motifs[index] = profileProbable(DNA[index], k, profile)
        if score(motifs) < score(bestMotifs):
            bestMotifs = motifs.copy()
            count +=1
        else:
            print(count)
            return bestMotifs

k = 8
t = 5
DNA = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
best = randomMotifSearch(DNA, k, t)
min = score(best)
for index in range(1000):
    print(index)
    a = randomMotifSearch(DNA, k, t)
    if score(a) < score(best):
        best = a
        min = score(a)
print(min)
for index in best:
    print(index)