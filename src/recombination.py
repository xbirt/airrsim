"""
recombination.py

This module simulates the V(D)J recombination process for adaptive immune receptors.
It includes functions for adding P and N nucleotides, performing exonuclease trimming,
and recombining V, D, and J segments to create clonotypes.

The module also includes a function for applying somatic hypermutation (SHM) to the
recombined sequences, which is relevant for B cell receptors.

Key functions:
- add_p_nucleotides: Adds P-nucleotides to segment ends
- add_n_nucleotides: Adds random N-nucleotides
- recombine_segments: Performs the V(D)J recombination process
- apply_shm: Applies somatic hypermutation to a sequence

This module is a crucial part of the adaptive immune receptor repertoire simulation,
providing a detailed and biologically relevant model of the recombination process.
"""

import random
import numpy as np
import re

SHM_HOTSPOTS = {
    'WRC': 5,    # W = A/T, R = A/G
    'GYW': 5,    # Y = C/T
    'WRCY': 4,
    'RGYW': 4,
    'WA': 1.5,
    'TW': 1.5,
    'SYC': 3,    # S = C/G
    'GRS': 3,
}

NUCLEOTIDE_MAP = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'W': '[AT]', 'R': '[AG]', 'Y': '[CT]', 'S': '[CG]',
    'H': '[ACT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
    'V': '[ACG]', 'D': '[AGT]', 'N': '[ACGT]'
}

# Pre-compile regex patterns for each hotspot motif
HOTSPOT_PATTERNS = {
    motif: re.compile(''.join(NUCLEOTIDE_MAP[base] for base in motif))
    for motif in SHM_HOTSPOTS.keys()
}

def create_mutation_likelihood_array(sequence):
    """
    Create an array of mutation likelihoods for each position in the sequence.
    
    Parameters:
    sequence (str): The DNA sequence to analyze.
    
    Returns:
    list: An array of mutation likelihoods for each position.
    """
    likelihood_array = [1.0] * len(sequence)  # Initialize with base likelihood of 1

    for motif, likelihood in SHM_HOTSPOTS.items():
        for match in HOTSPOT_PATTERNS[motif].finditer(sequence):
            start = match.start()
            if likelihood > likelihood_array[start]:
                likelihood_array[start] = likelihood  # Increase likelihood for the first position of the motif

    return likelihood_array

def perform_shm(sequence, mutation_rate=0.001):
    """
    Perform Somatic Hypermutation (SHM) on the given sequence, considering hotspots.
    Ensures each position is mutated at most once.
    
    Parameters:
    sequence (str): The DNA sequence to mutate.
    mutation_rate (float): The base mutation rate for non-hotspot positions.
    
    Returns:
    tuple: The mutated sequence and the number of mutations applied.
    """
    mutation_likelihoods = create_mutation_likelihood_array(sequence)
    expected_mutations = len(sequence) * mutation_rate
    
    # Calculate actual number of mutations based on Poisson distribution
    mutation_count = min(np.random.poisson(expected_mutations), len(sequence))
    
    sequence_list = list(sequence)
    available_positions = list(range(len(sequence)))
    
    for _ in range(mutation_count):
        # Choose a position based on weighted probabilities of available positions
        position = random.choices(available_positions, 
                                  weights=[mutation_likelihoods[i] for i in available_positions])[0]
        
        original_base = sequence_list[position]
        sequence_list[position] = random.choice([b for b in "ACGT" if b != original_base])
        
        available_positions.remove(position)

    return ''.join(sequence_list), mutation_count

def get_p_nucleotides(sequence, max_length=2):
    """
    Generate P-nucleotides based on the given sequence.
    
    Parameters:
    sequence (str): The DNA sequence to process.
    max_length (int): Maximum number of P-nucleotides to generate (usually 1-2).
    
    Returns:
    str: Generated P-nucleotides.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    p_length = random.randint(0, max_length)
    return ''.join(complement[base] for base in reversed(sequence[:p_length]))

def get_n_nucleotides(max_length=12):
    """
    Add N-nucleotides (random nucleotides).
    
    Parameters:
    max_length (int): Maximum number of N-nucleotides to add.
    
    Returns:
    str: String of random nucleotides.
    """
    return ''.join(random.choice("ATGC") for _ in range(random.randint(0, max_length)))

def recombine_segments(v, d, j, apply_shm=False, shm_rate=0.001):
    """
    Recombine V, (D), and J segments to create a clonotype.
    Optionally applies Somatic Hypermutation.
    
    Parameters:
    v, d, j (str): Coding sequences of V, D, and J segments.
    apply_shm (bool): Whether to apply Somatic Hypermutation.
    shm_rate (float): The mutation rate for SHM.
    
    Returns:
    str: Recombined sequence representing a clonotype, potentially with SHM.
    """
    max_trim_v_j = 15
    recombination_info = {}

    # Exonuclease trimming V segment
    v_trim_length = random.randint(0, min(len(v), max_trim_v_j))
    v_end = v[:-v_trim_length] if v_trim_length > 0 else v
    v_p_nucleotides = get_p_nucleotides(v_end)
    v_end += v_p_nucleotides
    recombination_info['v_region_len'] = len(v)
    recombination_info['v_region_trim_len'] = v_trim_length
    recombination_info['v_region_p_len'] = len(v_p_nucleotides)

    # Exonuclease trimming J segment
    j_trim_length = random.randint(0, min(len(j), max_trim_v_j))
    j_start = j[j_trim_length:] if j_trim_length > 0 else j
    j_p_nucleotides = get_p_nucleotides(j_start)
    j_start = j_p_nucleotides + j_start
    recombination_info['j_region_len'] = len(j)
    recombination_info['j_region_trim_len'] = j_trim_length
    recombination_info['j_region_p_len'] = len(j_p_nucleotides)

    if d is not None and d != '':
        # Exonuclease trimming D segment
        d_5_trim = random.randint(0, len(d) // 2)
        d_3_trim = random.randint(0, len(d) - d_5_trim - 1)
        d_mid = d[d_5_trim:len(d)-d_3_trim]
        d_p_nucleotides_5 = get_p_nucleotides(d_mid[::-1])
        d_p_nucleotides_3 = get_p_nucleotides(d_mid)
        d_mid = d_p_nucleotides_5[::-1] + d_mid + d_p_nucleotides_3

        # Add N-nucleotides and combine segments
        vd_junction = get_n_nucleotides(max_length=12)
        dj_junction = get_n_nucleotides(max_length=12)
        recombined_sequence = v_end + vd_junction + d_mid + dj_junction + j_start

        recombination_info['d_region_len'] = len(d)
        recombination_info['d_5_trim_len'] = d_5_trim
        recombination_info['d_3_trim_len'] = d_3_trim
        recombination_info['vd_junction_len'] = len(vd_junction)
        recombination_info['dj_junction_len'] = len(dj_junction)
    else:
        # If no D segment or empty D segment, just add N-nucleotides between V and J
        vj_junction = get_n_nucleotides(max_length=15)  # Slightly more for V-J junctions
        recombined_sequence = v_end + vj_junction + j_start
        recombination_info['vj_junction_len'] = len(vj_junction)

    # Apply SHM to the V region and V-D-J junction if specified
    shm_count = 0
    if apply_shm:
        v_and_junction = recombined_sequence[:-len(j_start)]
        mutated_v_and_junction, shm_count = perform_shm(v_and_junction, shm_rate)
        recombined_sequence = mutated_v_and_junction + j_start
    
    recombination_info['shm_count'] = shm_count

    return recombined_sequence, recombination_info