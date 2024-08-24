"""
read_simulation.py

This module handles the simulation of sequencing reads from clonotype sequences.
It provides functionality to generate reads of specified length from a FASTA file
containing clonotypes, with options to avoid reads entirely within the constant region.
"""

import os
import random
import re
import warnings
from src.utils import save_fasta, read_fasta_with_error_handling

CLONOTYPE_MAX_RETRIES = 3

def parse_sequence_id(seq_id):
    """
    Parse the sequence ID to extract J and C region lengths if possible.
    
    Args:
    seq_id (str): The sequence ID string.
    
    Returns:
    tuple: (j_length, c_length) or (None, None) if not found
    """
    parts = seq_id.split('_')
    j_info = next((part for part in parts if part.startswith('J')), None)
    c_info = next((part for part in parts if part.startswith('C')), None)
    
    if j_info and c_info:
        j_length = sum(int(x) for x in re.findall(r'[+-]?\d+', j_info))
        c_length = int(c_info[1:]) if c_info else 0
        return j_length, c_length
    return None, None

def generate_reads(input_file, output_file, read_length, read_count, with_constant_region=True):
    """
    Generate simulated reads from clonotype sequences.
    
    Args:
    input_file (str): Path to the input FASTA file containing clonotypes.
    output_file (str): Path to the output FASTA file for simulated reads.
    read_length (int): Length of each simulated read.
    read_count (int): Number of reads to generate.
    with_constant_region (bool): If True, include reads from the constant region. If False, try to avoid generating reads entirely within the constant region.
    
    Returns:
    None
    """
    # Used when avoiding the constant region
    minimum_overlap_left_of_j_region = 12

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    try:
        clonotypes = read_fasta_with_error_handling(input_file, "Input clonotype")
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {str(e)}")
        print("Read simulation aborted.")
        return
    
    # Check if we can avoid constant region
    sample_id = clonotypes[0][0]
    j_length, c_length = parse_sequence_id(sample_id)
    if not with_constant_region and (j_length is None or c_length is None):
        warnings.warn("Cannot avoid constant region due to incompatible sequence ID format. Proceeding with reads from all regions.")
        with_constant_region = True
    
    def read_generator():
        for _ in range(read_count):
            attempts = 0
            while True:
                try:
                    clonotype = random.choice(clonotypes)
                    seq_id, seq = clonotype
                    if len(seq) < read_length:
                        print(f"Warning: Sequence length ({len(seq)}) is smaller than read length ({read_length}) for sequence ID: {seq_id}. Trying to pick another clonotype.")
                        raise ValueError("Sequence too short")
                    
                    break

                except ValueError as e:
                    attempts += 1
                    if attempts == CLONOTYPE_MAX_RETRIES:
                        raise ValueError(f"Failed to find a suitable clonotype after 10 attempts: {str(e)}")
                    
            if not with_constant_region:
                j_length, c_length = parse_sequence_id(seq_id)
                max_start = len(seq) - c_length - j_length - minimum_overlap_left_of_j_region
                if max_start < 0:
                    raise ValueError(f"Sequence is too short to generate reads with the given parameters for sequence ID: {seq_id}")

                max_start = min(len(seq) - read_length, max_start)
                start_pos = random.randint(0, max_start)
            else:
                start_pos = random.randint(0, len(seq) - read_length)
            
            read_seq = seq[start_pos:start_pos + read_length]
            
            # Modify the sequence ID format
            parts = seq_id.split('_')
            clonotype_index = next((i for i, part in enumerate(parts) if part.startswith('clonotype')), None)
            
            if clonotype_index is not None and clonotype_index + 1 < len(parts):
                receptor_type = parts[0]
                clonotype_id = parts[clonotype_index + 1]  # Get the part after 'clonotype'
                new_id_start = f"{receptor_type}_c{clonotype_id}_R{start_pos}"
                new_id_rest = '_'.join(parts[1:clonotype_index] + parts[clonotype_index+2:])
                read_id = f"{new_id_start}_{new_id_rest}"
            else:
                read_id = f"{seq_id}_R{start_pos}"
            
            yield (read_id, read_seq)

    save_fasta(read_generator(), output_file)

    c_region_message = "including reads from all regions" if with_constant_region else "avoiding reads entirely within the constant region"
    print(f"Generated {read_count} reads of length {read_length}, {c_region_message}, and saved to {output_file}")