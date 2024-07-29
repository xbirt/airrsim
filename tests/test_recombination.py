import unittest
from src.recombination import (create_mutation_likelihood_array, perform_shm, 
                           get_p_nucleotides, get_n_nucleotides, recombine_segments)

class TestRecombination(unittest.TestCase):

    def test_create_mutation_likelihood_array(self):
        sequence = "AGTCTGCGCTTAA"
        likelihood_array = create_mutation_likelihood_array(sequence)
        self.assertEqual(len(likelihood_array), len(sequence))
        self.assertGreater(likelihood_array[4], 1.0)  # 'WRC' hotspot
        self.assertGreater(likelihood_array[7], 1.0)  # 'GYW' hotspot
        self.assertEqual(likelihood_array[0], 1.0)    # No hotspot

    def test_perform_shm(self):
        sequence = "AGTCATGCATGCATGC"
        mutated_seq, mutation_count = perform_shm(sequence, mutation_rate=0.1)
        self.assertEqual(len(mutated_seq), len(sequence))
        self.assertNotEqual(mutated_seq, sequence)
        self.assertGreater(mutation_count, 0)

    def test_get_p_nucleotides(self):
        sequence = "AGTC"
        p_nucleotides = get_p_nucleotides(sequence, max_length=2)
        self.assertLessEqual(len(p_nucleotides), 2)
        for i, base in enumerate(p_nucleotides):
            self.assertIn(base, "ACGT")
            self.assertNotEqual(base, sequence[i])

    def test_get_n_nucleotides(self):
        n_nucleotides = get_n_nucleotides(max_length=12)
        self.assertLessEqual(len(n_nucleotides), 12)
        self.assertTrue(all(base in "ACGT" for base in n_nucleotides))

    def test_recombine_segments_with_d(self):
        v = "CGTCGTCGTCGT"
        d = "ATATATATATAT"
        j = "GCAGCAGCAGCA"
        recombined, info = recombine_segments(v, d, j)
        self.assertIn('v_region_trim_len', info)
        self.assertIn('d_5_trim_len', info)
        self.assertIn('d_3_trim_len', info)
        self.assertIn('j_region_trim_len', info)
        self.assertIn('vd_junction_len', info)
        self.assertIn('dj_junction_len', info)

    def test_recombine_segments_without_d(self):
        v = "CGTCGTCGTCGT"
        j = "GCAGCAGCAGCA"
        recombined, info = recombine_segments(v, None, j)
        self.assertIn('v_region_trim_len', info)
        self.assertIn('j_region_trim_len', info)
        self.assertIn('vj_junction_len', info)
        self.assertNotIn('d_5_trim_len', info)

    def test_recombine_segments_with_shm(self):
        v = "CGTCGTCGTCGT"
        d = "ATATATATATAT"
        j = "GCAGCAGCAGCA"
        recombined, info = recombine_segments(v, d, j, apply_shm=True, shm_rate=1)
        self.assertIn('shm_count', info)
        self.assertGreater(info['shm_count'], 0)

    def test_edge_cases(self):
        # Test with very short segments
        v, d, j = "A", "T", "G"
        recombined, info = recombine_segments(v, d, j)
        self.assertGreater(len(recombined), 3)

        # Test with empty segments
        v, d, j = "", "", ""
        recombined, info = recombine_segments(v, d, j)
        self.assertGreater(len(recombined), 0)

    def test_consistency(self):
        v = "CGTCGTCGTCGT"
        d = "ATATATATATAT"
        j = "GCAGCAGCAGCA"
        recombined1, _ = recombine_segments(v, d, j, apply_shm=False)
        recombined2, _ = recombine_segments(v, d, j, apply_shm=False)
        self.assertNotEqual(recombined1, recombined2)  # Random elements should make these different

if __name__ == '__main__':
    unittest.main()