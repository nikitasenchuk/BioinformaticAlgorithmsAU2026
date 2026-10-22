"""Запуск решения BA2C или BA2E на текстовом файле Rosalind."""
import argparse
from pathlib import Path
from motifs import NUCLEOTIDES, greedy_motif_search, profile_most_probable_kmer


def solve(problem: str, data: str) -> str:
    lines = [line.strip() for line in data.splitlines() if line.strip()]
    if problem == 'ba2c':
        if len(lines) != 6:
            raise ValueError('BA2C: ожидаются Text, k и четыре строки профиля.')
        text, k = lines[0], int(lines[1])
        # В файлах Rosalind строки матрицы идут в порядке A, C, G, T.
        profile = {base: list(map(float, row.split()))
                   for base, row in zip(NUCLEOTIDES, lines[2:])}
        return profile_most_probable_kmer(text, k, profile) + '\n'
    if problem == 'ba2e':
        k, t = map(int, lines[0].split())
        dna = ' '.join(lines[1:]).split()  # Наборы могут разделять DNA пробелами или строками.
        if len(dna) != t:
            raise ValueError(f'Ожидалось {t} последовательностей, получено {len(dna)}.')
        return '\n'.join(greedy_motif_search(dna, k)) + '\n'
    raise ValueError('Неизвестная задача: ' + problem)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('problem', choices=['ba2c', 'ba2e'])
    parser.add_argument('input', type=Path)
    parser.add_argument('-o', '--output', type=Path)
    args = parser.parse_args()
    answer = solve(args.problem, args.input.read_text(encoding='utf-8-sig'))
    if args.output:
        args.output.write_text(answer, encoding='utf-8')
    else:
        print(answer, end='')


if __name__ == '__main__':
    main()
