"""
utils.py

This module provides utility functions for working with receptor sequences and FASTA files.
It includes functionality for:
- Determining if a receptor type has a D segment
- Reading FASTA files with error handling
- Saving sequences to FASTA format

The module relies on the Biopython library for sequence handling and provides
robust error checking and file I/O operations.

Functions:
- receptor_has_d_segment: Check if a receptor type has a D segment
- read_fasta_with_error_handling: Read a FASTA file with comprehensive error handling
- read_fasta: Read and yield sequences from a FASTA file
- save_fasta: Save sequences to a FASTA file with specified line length
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def receptor_has_d_segment(receptor_type):
    """
    Determine if a receptor type has a D segment.

    Args:
    receptor_type (str): The type of receptor (e.g., 'TRA', 'TRB', 'IGH').

    Returns:
    bool: True if the receptor type has a D segment, False otherwise.
    """
    # Define which receptor types have D segments
    RECEPTORS_WITH_D = ['IGH', 'TRB', 'TRD']

    return receptor_type in RECEPTORS_WITH_D

def read_fasta_with_error_handling(file_path, segment_type=""):
    """
    Read a FASTA file with error handling.

    Parameters:
    file_path (str): Path to the FASTA file.
    segment_type (str, optional): Type of segment (e.g., 'V', 'D', 'J', 'C') for error messaging.

    Returns:
    list: List of tuples containing (header, sequence) from the FASTA file.

    Raises:
    FileNotFoundError: If the file does not exist.
    ValueError: If the file is empty or not in the correct FASTA format.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{segment_type} file not found: {file_path}")

    try:
        segments = list(read_fasta(file_path))
        if not segments:
            raise ValueError(f"No sequences found in the {segment_type} file: {file_path}")
        return segments
    except Exception as e:
        raise ValueError(f"Error reading the {segment_type} file. Please ensure it is a valid FASTA format: {str(e)}")

def read_fasta(file_path):
    """
    Read sequences from a FASTA file and yield them as tuples.

    Parameters:
    file_path (str): The path to the FASTA file to be read.

    Yields:
    tuple: A tuple containing two items:
           1. The sequence identifier (str)
           2. The sequence (str)

    Example:
    >>> for seq_id, sequence in read_fasta("path/to/your/file.fasta"):
    ...     print(seq_id, len(sequence))
    """
    if os.path.exists(file_path):
        for record in SeqIO.parse(file_path, "fasta"):
            yield (record.id, str(record.seq).upper())
    else:
        print(f"Warning: File {file_path} does not exist.")
        return iter([])  # Return an empty iterator

def save_fasta(sequences, file_path, line_length=60):
    """
    Save a sequence of (id, sequence) tuples to a FASTA file with specified line length.

    Parameters:
    sequences (iterable of tuple): An iterable where each element is a tuple containing:
                                   1. The sequence identifier (str)
                                   2. The sequence (str)
    file_path (str): The path where the FASTA file will be saved.
    line_length (int): The number of characters per line for the sequence. Default is 60.

    Example:
    >>> save_fasta(generate_sequences(), "path/to/output.fasta", line_length=80)
    """
    def sequence_generator():
        for seq_id, seq in sequences:
            yield SeqRecord(Seq(seq), id=seq_id, description="")

    with open(file_path, "w") as output_handle:
        writer = SeqIO.FastaIO.FastaWriter(output_handle, wrap=line_length)
        writer.write_file(sequence_generator())