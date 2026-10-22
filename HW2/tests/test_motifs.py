"""Проверки примеров Rosalind, граничных случаев и правил разрешения равенств."""
import itertools
from pathlib import Path
import random
import unittest
from motifs import (score, profile_with_pseudocounts, profile_probability,
                    profile_most_probable_kmer, greedy_motif_search)
from run_rosalind import solve

DATA = Path(__file__).parent / 'data'

class MotifTests(unittest.TestCase):
    def test_rosalind_samples(self):
        for problem in ('ba2c', 'ba2e'):
            with self.subTest(problem=problem):
                self.assertEqual(solve(problem, (DATA / f'{problem}_sample.txt').read_text()),
                                 (DATA / f'{problem}_expected.txt').read_text())

    def test_score_sample(self):
        self.assertEqual(score(['TTC', 'ATC', 'TTC', 'ATC', 'TTC']), 2)

    def test_score_against_exhaustive_consensus(self):
        rng = random.Random(2026)
        for _ in range(80):
            motifs = [''.join(rng.choices('ACGT', k=3)) for _ in range(5)]
            # Независимая проверка: перебираем все возможные consensus длины 3.
            expected = min(sum(a != b for motif in motifs for a, b in zip(motif, consensus))
                           for consensus in itertools.product('ACGT', repeat=3))
            self.assertEqual(score(motifs), expected)

    def test_pseudocount_values_and_normalization(self):
        profile = profile_with_pseudocounts(['AA', 'AC'])
        self.assertEqual(profile['A'], [3 / 6, 2 / 6])
        self.assertEqual(profile['C'], [1 / 6, 2 / 6])
        for col in range(2):
            self.assertAlmostEqual(sum(row[col] for row in profile.values()), 1.0)
            self.assertTrue(all(row[col] > 0 for row in profile.values()))

    def test_probability(self):
        profile = {'A': [0.2, 0.3], 'C': [0.3, 0.2], 'G': [0.4, 0.4], 'T': [0.1, 0.1]}
        self.assertAlmostEqual(profile_probability('AC', profile), 0.04)

    def test_first_window_on_equal_probabilities(self):
        uniform = {base: [0.25] * 2 for base in 'ACGT'}
        self.assertEqual(profile_most_probable_kmer('TGAC', 2, uniform), 'TG')

    def test_first_window_when_all_candidates_have_probability_zero(self):
        profile = {'A': [1.0], 'C': [0.0], 'G': [0.0], 'T': [0.0]}
        self.assertEqual(profile_most_probable_kmer('TGCC', 1, profile), 'T')

    def test_last_window_is_considered(self):
        profile = {base: [float(base == 'G')] for base in 'ACGT'}
        self.assertEqual(profile_most_probable_kmer('AAAG', 1, profile), 'G')

    def test_greedy_preserves_initial_set_on_tie(self):
        self.assertEqual(greedy_motif_search(['TA', 'AT'], 1), ['T', 'T'])
        self.assertEqual(greedy_motif_search(['TA', 'CG'], 1), ['T', 'C'])

    def test_single_sequence(self):
        self.assertEqual(greedy_motif_search(['GACT'], 2), ['GA'])

    def test_k_equals_sequence_length(self):
        self.assertEqual(greedy_motif_search(['ACG', 'ATG'], 3), ['ACG', 'ATG'])

    def test_input_is_not_modified(self):
        dna = ['ACGT', 'TGCA']
        before = dna.copy()
        greedy_motif_search(dna, 2)
        self.assertEqual(dna, before)

    def test_invalid_inputs(self):
        self.assertEqual(score([]), 0)
        for motifs in ([], [''], ['AA', 'A'], ['AN']):
            with self.subTest(motifs=motifs), self.assertRaises(ValueError):
                profile_with_pseudocounts(motifs)
        for dna, k in (([], 1), (['A'], 2), (['A'], 0), (['AN'], 1)):
            with self.subTest(dna=dna, k=k), self.assertRaises(ValueError):
                greedy_motif_search(dna, k)

    def test_cli_rejects_wrong_sequence_count(self):
        with self.assertRaises(ValueError):
            solve('ba2e', '2 3\nACGT\nTGCA\n')

if __name__ == '__main__':
    unittest.main()
