"""
simulation.py

This module contains functions for simulating adaptive immune receptor repertoires.
It provides functionality to generate clonotypes for various receptor types
(e.g., immunoglobulins and T cell receptors) based on input gene segment files.

Key components:
1. Clonotype generation: Randomly selects V, (D), and J segments and recombines them.
2. Constant region handling: Adds appropriate constant regions to the recombined sequences.
3. Somatic Hypermutation (SHM): Applies SHM to B cell receptor sequences if specified.
4. FASTA file I/O: Reads input gene segments and writes generated repertoires to FASTA files.

Main functions:
- generate_clonotypes: Generates individual clonotypes for a given receptor type.
- simulate_receptor_repertoires: Simulates entire repertoires and saves them to files.

Helper functions:
- extract_segment_id: Extracts segment IDs from FASTA headers.
- load_constant_regions: Loads constant region sequences from files.
- get_constant_region: Selects appropriate constant regions for clonotypes.

This module is designed to work with IMGT-formatted gene segment files and can
generate large-scale repertoires efficiently through streaming processing.

Usage example:
    simulate_receptor_repertoires('IGH', 'output_igh_repertoire.fasta',
                                  'IGHV.fasta', 'IGHJ.fasta', 'IGHD.fasta', 'IGHC.fasta',
                                  num_clonotypes=1000000, shm_rate=0.002)

Note: This simulation provides a simplified model of the biological VDJ recombination
process and is intended for generating synthetic datasets for immunoinformatics research.
"""

import random
import re
from src.utils import read_fasta, save_fasta, read_fasta_with_error_handling
from src.recombination import recombine_segments

def extract_segment_id(header):
    """
    Extract the segment ID from the FASTA header.
    
    Parameters:
    header (str): The FASTA header line.
    
    Returns:
    str: The extracted segment ID.
    """
    match = re.search(r'[^|]*\|([^|]+)\|', header)
    if match:
        return match.group(1)
    return "unknown"

def load_constant_regions(c_file):
    """
    Load constant regions from a given file.
    
    Parameters:
    c_file (str): The path to the constant region file.
    
    Returns:
    list: A list of tuples containing (constant_region_id, constant_region_sequence).
    """
    try:
        constant_regions = list(read_fasta(c_file))
        if not constant_regions:
            raise ValueError(f"No constant region sequences found in file: {c_file}")
        return constant_regions
    except Exception as e:
        raise ValueError(f"Error reading the constant region file. Please ensure it is a valid FASTA format: {str(e)}")

def get_constant_region(receptor_type, constant_regions):
    """
    Select an appropriate constant region for the given receptor type.
    
    Parameters:
    receptor_type (str): The type of receptor (e.g., 'IGH', 'IGL', 'IGK', 'TRA', 'TRB').
    constant_regions (list): Pre-loaded list of constant regions.
    
    Returns:
    tuple: A tuple containing the selected constant region sequence and its ID.
    
    Raises:
    ValueError: If the receptor type is unknown or if no matching constant region is found.
    """
    if receptor_type == 'IGH':
        # For IGH, implement a simplified class switching model
        isotype_weights = {'IGHM': 0.4, 'IGHD': 0.1, 'IGHG': 0.3, 'IGHA': 0.15, 'IGHE': 0.05}
        selected_isotype = random.choices(list(isotype_weights.keys()), 
                                          weights=list(isotype_weights.values()))[0]
        
        # Find the constant region matching the selected isotype
        for const_id, const_seq in constant_regions:
            if extract_segment_id(const_id).startswith(selected_isotype):
                return const_seq, extract_segment_id(const_id)
        
        raise ValueError(f"No matching constant region found for {selected_isotype}")
    
    elif receptor_type in ['IGL', 'IGK', 'TRA', 'TRB', 'TRD', 'TRG']:
        # For other receptor types, randomly choose between available constant regions
        matching_regions = [cr for cr in constant_regions if extract_segment_id(cr[0]).startswith(receptor_type)]
        if not matching_regions:
            raise ValueError(f"No matching constant region found for {receptor_type}")
        const_id, const_seq = random.choice(matching_regions)
        return const_seq, extract_segment_id(const_id)
    
    else:
        raise ValueError(f"Unknown receptor type: {receptor_type}")

def generate_clonotypes(receptor_type, num_clonotypes, v_file, j_file, d_file=None,
                         c_file=None, apply_shm=False, shm_rate=0.001,
                         append_constant_region=True, preserve_alignment=True):
    """
    Generate clonotypes for a given receptor type.

    Parameters:
    receptor_type (str): The type of receptor to simulate (e.g., 'IGH', 'TRA', 'TRB').
    num_clonotypes (int): The number of clonotypes to generate.
    v_file (str): Path to the V gene segment file.
    j_file (str): Path to the J gene segment file.
    d_file (str, optional): Path to the D gene segment file (if applicable).
    c_file (str, optional): Path to the constant region file (if applicable).
    apply_shm (bool): Whether to apply Somatic Hypermutation (only for B cell receptors).
    shm_rate (float): The mutation rate for SHM.
    include_constant (bool): Whether to include the constant region in the output sequences.
    preserve_alignment (bool): Whether to preserve the codon alignment when performing the recombination.

    Yields:
    tuple: A tuple containing the clonotype ID and sequence.
    """
    try:
        # Read V and J segments
        v_segments = read_fasta_with_error_handling(v_file, 'V')
        j_segments = read_fasta_with_error_handling(j_file, 'J')
        
        # Read D segments if applicable
        d_segments = read_fasta_with_error_handling(d_file, 'D') if d_file else None

        # Load constant regions if needed
        constant_regions = None
        if append_constant_region and c_file:
            constant_regions = load_constant_regions(c_file)
        
        for i in range(num_clonotypes):
            # Randomly select V and J segments
            v_segment = random.choice(v_segments)
            j_segment = random.choice(j_segments)
            
            # Randomly select D segment if applicable
            d_segment = random.choice(d_segments) if d_segments else None

            # Extract segment IDs
            v_id = extract_segment_id(v_segment[0])
            j_id = extract_segment_id(j_segment[0])
            d_id = extract_segment_id(d_segment[0]) if d_segment else None
            
            # Recombine segments to create a clonotype
            variable_region, recombination_info = recombine_segments(
                v_segment[1], d_segment[1] if d_segment else None, j_segment[1],
                apply_shm=apply_shm, shm_rate=shm_rate, preserve_alignment=preserve_alignment
            )
            
            # Add constant region if specified
            c_id = None
            if append_constant_region and constant_regions:
                constant_region, c_id = get_constant_region(receptor_type, constant_regions)
                clonotype_seq = variable_region + constant_region
                c_region_len = len(constant_region)
            else:
                clonotype_seq = variable_region

            # Resolve ambiguous nucleotides, if any
            clonotype_seq = resolve_ambiguous_nucleotides(clonotype_seq)

            # Create an informative clonotype ID
            clonotype_id_parts = [receptor_type, f"clonotype_{i+1}", v_id]
            if d_id:
                clonotype_id_parts.append(d_id)
            clonotype_id_parts.append(j_id)
            if c_id:
                clonotype_id_parts.append(c_id)
            
            # Add recombination information to clonotype ID
            v_info = f"V{recombination_info['v_region_len']}-{recombination_info['v_region_trim_len']}+{recombination_info['v_region_p_len']}"
            
            vd_info = ''
            vj_info = ''
            d_info = ''
            if d_id:
                vd_info = f"+{recombination_info['vd_junction_len']}"
                d_info = f"_D{recombination_info['d_region_len']}-{recombination_info['d_5_trim_len']}-{recombination_info['d_3_trim_len']}"
                j_info = f"J{recombination_info['dj_junction_len']}+{recombination_info['j_region_p_len']}+{recombination_info['j_region_len']}-{recombination_info['j_region_trim_len']}"
            else:
                vj_info = f"+{recombination_info['vj_junction_len']}"
                j_info = f"J{recombination_info['j_region_p_len']}+{recombination_info['j_region_len']}-{recombination_info['j_region_trim_len']}"
            
            clonotype_id_parts.append(f"{v_info}{vd_info}{vj_info}{d_info}")
            clonotype_id_parts.append(f"{j_info}")
            
            if append_constant_region and constant_regions:
                clonotype_id_parts.append(f"C{c_region_len}")
            
            if recombination_info['shm_count'] > 0:
                clonotype_id_parts.append(f"SHM{recombination_info['shm_count']}")
            
            clonotype_id = "_".join(clonotype_id_parts)

            # Create a tuple with clonotype ID and sequence
            yield (clonotype_id, clonotype_seq)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {str(e)}")
        print("Clonotype generation aborted.")
        return

def resolve_ambiguous_nucleotides(sequence):
    """
    Replace ambiguous nucleotides with specific nucleotides based on IUPAC notation.
    """
    iupac_map = {
        'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
        'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
        'N': 'ACGT'
    }
    
    resolved_sequence = ''
    for nucleotide in sequence.upper():
        if nucleotide in 'ACGT':
            resolved_sequence += nucleotide
        elif nucleotide in iupac_map:
            resolved_sequence += random.choice(iupac_map[nucleotide])
        else:
            raise ValueError(f"Unexpected character in sequence: {nucleotide}")
    
    return resolved_sequence

def simulate_receptor_repertoires(receptor_type, output_file, v_file, j_file, d_file=None, c_file=None, 
                                  num_clonotypes=1000, shm_rate=0.001, append_constant_region=True, apply_shm=True,
                                  preserve_alignment=True):
    """
    Simulate a receptor repertoire and save it to a FASTA file.

    This function generates clonotypes and saves them directly to a FASTA file,
    processing them in a streaming fashion to handle large numbers efficiently.
    It applies SHM to B cell receptors (IG types) if apply_shm is True and shm_rate > 0.

    Parameters:
    receptor_type (str): The type of receptor to simulate (e.g., 'IGH', 'TRA', 'TRB').
    output_file (str): The path to save the output FASTA file.
    v_file (str): Path to the V gene segment file.
    j_file (str): Path to the J gene segment file.
    d_file (str, optional): Path to the D gene segment file (if applicable).
    c_file (str, optional): Path to the constant region file (if applicable).
    num_clonotypes (int): The number of clonotypes to simulate. Default is 1000.
    shm_rate (float): The mutation rate for SHM (only applies to IG types). Default is 0.001 (0.1%).
    include_constant (bool): Whether to include the constant region in the output sequences.
    apply_shm (bool): Whether to apply somatic hypermutation. Default is True.
    preserve_alignment (bool): Whether to preserve the codon alignment when performing the recombination.

    Example:
    >>> simulate_receptor_repertoires('IGH', 'output_igh_repertoire.fasta', 'IGHV.fasta', 'IGHJ.fasta', 'IGHD.fasta', 'IGHC.fasta', num_clonotypes=1000000, shm_rate=0.002)
    >>> simulate_receptor_repertoires('TRA', 'output_tra_repertoire.fasta', 'TRAV.fasta', 'TRAJ.fasta', c_file='TRAC.fasta', num_clonotypes=1000000)
    >>> simulate_receptor_repertoires('IGH', 'output_igh_no_shm.fasta', 'IGHV.fasta', 'IGHJ.fasta', 'IGHD.fasta', num_clonotypes=1000000, shm_rate=0, apply_shm=False)
    """
    apply_shm_final = apply_shm and receptor_type.startswith('IG') and shm_rate > 0
    clonotypes = generate_clonotypes(receptor_type, num_clonotypes, v_file, j_file, d_file, c_file,
                                     apply_shm=apply_shm_final, shm_rate=shm_rate,
                                     append_constant_region=append_constant_region,
                                     preserve_alignment=preserve_alignment)
    save_fasta(clonotypes, output_file)
    
    print(f"Generated {num_clonotypes} {receptor_type} clonotypes" + 
          (f" with SHM rate {shm_rate}" if apply_shm_final else " without SHM") +
          (" including constant region" if append_constant_region else "") +
          f" and saved to {output_file}")