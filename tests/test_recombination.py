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

    def run_multiple_trials(self, test_function, num_trials=100):
        for _ in range(num_trials):
            test_function()

    def test_recombine_segments_alignment_preservation_with_d(self):
        def single_trial():
            v = "CGTCGTCGTCGT"  # 12 bp
            d = "ATATATATATAT"  # 12 bp
            j = "GCAGCAGCAGCA"  # 12 bp
            recombined, info = recombine_segments(v, d, j, preserve_alignment=True)
            
            # Check if the total length is divisible by 3
            self.assertEqual(len(recombined) % 3, 0, f"Length not divisible by 3: {recombined}")
            
            # Check if V and J are intact (allowing for trimming)
            self.assertTrue(recombined.startswith(v[:len(v)-info['v_region_trim_len']]),
                            f"V region mismatch: {recombined}")
            self.assertTrue(recombined.endswith(j[info['j_region_trim_len']:]),
                            f"J region mismatch: {recombined}")
            
            # Check if the junction lengths are recorded
            self.assertIn('vd_junction_len', info)
            self.assertIn('dj_junction_len', info)

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_alignment_preservation_without_d(self):
        def single_trial():
            v = "CGTCGTCGTCGT"  # 12 bp
            j = "GCAGCAGCAGCA"  # 12 bp
            recombined, info = recombine_segments(v, None, j, preserve_alignment=True)
            
            # Check if the total length is divisible by 3
            self.assertEqual(len(recombined) % 3, 0, f"Length not divisible by 3: {recombined}")
            
            # Check if V and J are intact (allowing for trimming)
            self.assertTrue(recombined.startswith(v[:len(v)-info['v_region_trim_len']]),
                            f"V region mismatch: {recombined}")
            self.assertTrue(recombined.endswith(j[info['j_region_trim_len']:]),
                            f"J region mismatch: {recombined}")
            
            # Check if the junction length is recorded
            self.assertIn('vj_junction_len', info)

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_without_alignment_preservation(self):
        def single_trial():
            v = "CGTCGTCGTCGT"  # 12 bp
            d = "ATATATATATAT"  # 12 bp
            j = "GCAGCAGCAGCA"  # 12 bp
            recombined, info = recombine_segments(v, d, j, preserve_alignment=False)
            
            # We just check if the function runs without errors and produces a sequence
            self.assertGreater(len(recombined), 0, f"Empty recombined sequence")

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_alignment_preservation_edge_cases(self):
        def single_trial():
            # Test with very short segments
            v, d, j = "A", "T", "G"
            recombined, info = recombine_segments(v, d, j, preserve_alignment=True)
            
            # Check if J segment starts at a position divisible by 3
            j_start = len(recombined) - (len(j) - info['j_region_trim_len'])
            self.assertEqual(j_start % 3, 0, f"J segment not aligned: {recombined}")
            
            # Test with empty segments
            v, d, j = "", "", ""
            recombined, info = recombine_segments(v, d, j, preserve_alignment=True)
            self.assertEqual(len(recombined) % 3, 0, f"Length not divisible by 3 for empty segments: {recombined}")

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_alignment_consistency(self):
        def single_trial():
            v = "CGTCGTCGTCGT"
            d = "ATATATATATAT"
            j = "GCAGCAGCAGCA"
            recombined1, _ = recombine_segments(v, d, j, preserve_alignment=True)
            recombined2, _ = recombine_segments(v, d, j, preserve_alignment=True)
            
            # Both should be divisible by 3, but likely different due to random elements
            self.assertEqual(len(recombined1) % 3, 0, f"Length of recombined1 not divisible by 3: {recombined1}")
            self.assertEqual(len(recombined2) % 3, 0, f"Length of recombined2 not divisible by 3: {recombined2}")
            # Note: We don't check if they're different, as there's a small chance they could be the same

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_j_frame_preservation(self):
        def single_trial():
            v = "CGTCGTCGTCGT"  # 12 bp
            d = "ATATATATATAT"  # 12 bp
            j = "GCAGCAGCAGCA"  # 12 bp
            
            recombined, info = recombine_segments(v, d, j, preserve_alignment=True)
            
            # Calculate the start of the J segment, excluding P-nucleotides
            j_start = len(recombined) - (len(j) - info['j_region_trim_len'])
            recombined_j = recombined[j_start:]
            
            # Find the start of the first full codon in the remaining J segment
            remaining_j = j[info['j_region_trim_len']:]
            first_full_codon_start = (3 - (info['j_region_trim_len'] % 3)) % 3
            
            # Calculate the position of the first full J codon in the recombined sequence
            first_full_j_codon_pos = j_start + first_full_codon_start
            
            # Check if the first full codon of J is aligned (starts at a position divisible by 3)
            self.assertEqual(first_full_j_codon_pos % 3, 0, 
                            f"First full codon of J not aligned. "
                            f"Recombined sequence: {recombined}\n"
                            f"J start: {j_start}, First full J codon position: {first_full_j_codon_pos}")
            
            # Check if the trimmed J segment is at the end
            self.assertTrue(recombined.endswith(remaining_j), 
                            f"Trimmed J segment not at the end. Recombined sequence: {recombined}")
            
            # Check if the trimmed V segment is at the start
            trimmed_v = v[:len(v)-info['v_region_trim_len']]
            self.assertTrue(recombined.startswith(trimmed_v),
                            f"Trimmed V segment not at the start. Recombined sequence: {recombined}")

        self.run_multiple_trials(single_trial)

    def test_recombine_segments_j_frame_preservation_no_d(self):
        def single_trial():
            v = "CGTCGTCGTCGT"  # 12 bp
            j = "GCAGCAGCAGCA"  # 12 bp
            
            recombined, info = recombine_segments(v, None, j, preserve_alignment=True)
            
            # Calculate the start of the J segment, excluding P-nucleotides
            j_start = len(recombined) - (len(j) - info['j_region_trim_len'])
            recombined_j = recombined[j_start:]
            
            # Find the start of the first full codon in the remaining J segment
            remaining_j = j[info['j_region_trim_len']:]
            first_full_codon_start = (3 - (info['j_region_trim_len'] % 3)) % 3
            
            # Calculate the position of the first full J codon in the recombined sequence
            first_full_j_codon_pos = j_start + first_full_codon_start
            
            # Check if the first full codon of J is aligned (starts at a position divisible by 3)
            self.assertEqual(first_full_j_codon_pos % 3, 0, 
                            f"First full codon of J not aligned. "
                            f"Recombined sequence: {recombined}\n"
                            f"J start: {j_start}, First full J codon position: {first_full_j_codon_pos}")
            
            # Check if the trimmed J segment is at the end
            self.assertTrue(recombined.endswith(remaining_j), 
                            f"Trimmed J segment not at the end. Recombined sequence: {recombined}")
            
            # Check if the trimmed V segment is at the start
            trimmed_v = v[:len(v)-info['v_region_trim_len']]
            self.assertTrue(recombined.startswith(trimmed_v),
                            f"Trimmed V segment not at the start. Recombined sequence: {recombined}")

            # Check if VJ junction length is recorded
            self.assertIn('vj_junction_len', info,
                        f"VJ junction length not recorded in info: {info}")

        self.run_multiple_trials(single_trial)

if __name__ == '__main__':
    unittest.main()