from math import prod

NUCLEOTIDES = 'ACGT'


def _validate_motifs(motifs: list[str]) -> int:
    if not motifs or not motifs[0]:
        raise ValueError('Нужен хотя бы один непустой мотив.')
    k = len(motifs[0])
    if any(len(motif) != k for motif in motifs):
        raise ValueError('Мотивы должны иметь одинаковую длину.')
    if any(base not in NUCLEOTIDES for motif in motifs for base in motif):
        raise ValueError('Допустимы только нуклеотиды A, C, G, T.')
    return k


def score(motifs: list[str]) -> int:
    if not motifs:
        return 0
    k = _validate_motifs(motifs)
    consensus = []
    for column in range(k):
        counts = {base: sum(motif[column] == base for motif in motifs)
                  for base in NUCLEOTIDES}
        consensus.append(max(NUCLEOTIDES, key=counts.get))
    return sum(base != consensus[column]
               for motif in motifs for column, base in enumerate(motif))


def profile_with_pseudocounts(motifs: list[str]) -> dict[str, list[float]]:
    k = _validate_motifs(motifs)
    counts = {base: [1] * k for base in NUCLEOTIDES}
    for motif in motifs:
        for column, base in enumerate(motif):
            counts[base][column] += 1
    denominator = len(motifs) + 4
    return {base: [count / denominator for count in counts[base]]
            for base in NUCLEOTIDES}

def profile_probability(pattern: str, profile: dict[str, list[float]]) -> float:
    return prod(profile[base][column] for column, base in enumerate(pattern))


def profile_most_probable_kmer(
    text: str, k: int, profile: dict[str, list[float]],
) -> str:
    if k < 1 or k > len(text):
        raise ValueError('k должно быть от 1 до длины последовательности.')
    if any(base not in NUCLEOTIDES for base in text):
        raise ValueError('Допустимы только нуклеотиды A, C, G, T.')
    if set(profile) != set(NUCLEOTIDES) or any(len(row) != k for row in profile.values()):
        raise ValueError('Профиль должен содержать четыре строки A/C/G/T длины k.')
    best_pattern = text[:k]
    best_probability = profile_probability(best_pattern, profile)
    for start in range(1, len(text) - k + 1):
        pattern = text[start:start + k]
        probability = profile_probability(pattern, profile)
        if probability > best_probability:
            best_pattern, best_probability = pattern, probability
    return best_pattern