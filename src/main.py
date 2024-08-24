"""
main.py

This is the main entry point for the Adaptive Immune Receptor Repertoire (AIRR) Simulation Tool.
It provides a command-line interface for various actions including generating clonotypes,
simulating sequencing reads, and downloading reference data.

The script uses argparse to handle command-line arguments and calls appropriate functions
based on the specified action.
"""

import argparse
import os
from sys import exit
from src.clonotype_simulation import simulate_receptor_repertoires
from src.download_data import download_data
from src.read_simulation import generate_reads
from src.utils import receptor_has_d_segment

def check_file_exists(file_path, file_type):
    """
    Check if a file exists and raise an error if it doesn't.

    Args:
    file_path (str): Path to the file to check.
    file_type (str): Type of the file (e.g., 'V', 'D', 'J', 'C') for error messaging.

    Raises:
    FileNotFoundError: If the file does not exist.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Error: {file_type} file not found: {file_path}")

def case_insensitive_group(value):
    """
    Convert input to uppercase, but keep 'all' as is.
    """
    return value.upper() if value.lower() != 'all' else 'all'

def parse_args():
    """
    Parse command-line arguments for the AIRR Simulation Tool.

    Returns:
    argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Adaptive Immune Receptor Repertoire (AIRR) Simulation Tool")
    subparsers = parser.add_subparsers(dest='action', required=True, help='Action to perform')

    # Common arguments for generateClonotypes and downloadData
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('group', type=case_insensitive_group, choices=['TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGK', 'IGL', 'all'], help='Receptor group')
    common_parser.add_argument('--species', choices=['human', 'mouse', 'all'], default='human', help='Species (default: human)')

    # generateClonotypes subparser
    parser_generate = subparsers.add_parser('generateClonotypes', parents=[common_parser], help='Generate clonotypes')
    parser_generate.add_argument('--output', type=str, help='Output file (default: clonotypes[group].fasta)')
    parser_generate.add_argument('--count', type=int, default=1000, help='Number of clonotypes to generate')
    parser_generate.add_argument('--shm-rate', type=float, default=0.001, help='Somatic hypermutation rate')
    parser_generate.add_argument('--with-shm', action='store_true', default=True, help='Apply somatic hypermutations (default: True)')
    parser_generate.add_argument('--no-shm', action='store_false', dest='with_shm', help='Do not apply somatic hypermutations')
    constant_region_group = parser_generate.add_mutually_exclusive_group()
    constant_region_group.add_argument('--with-constant-region', '--include-constant-region', '--with-c-region', '--include-c-region', dest='with_constant_region', action='store_true', default=True, help='Include constant region in the output sequences (default)')
    constant_region_group.add_argument('--without-constant-region', '--exclude-constant-region', '--without-c-region', '--exclude-c-region', '--no-c-region', dest='with_constant_region', action='store_false', help='Exclude constant region from the output sequences')
    parser_generate.add_argument('--with-alignment-preservation', action='store_true', default=True, help='Preserve codon alignment in recombination (default: True)')
    parser_generate.add_argument('--without-alignment-preservation', action='store_false', dest='with_alignment_preservation', help='Do not preserve codon alignment in recombination')

    parser_generate.add_argument('--v_file', type=str, help='Path to V gene segment file')
    parser_generate.add_argument('--d_file', type=str, help='Path to D gene segment file')
    parser_generate.add_argument('--j_file', type=str, help='Path to J gene segment file')
    parser_generate.add_argument('--c_file', type=str, help='Path to C gene segment file')
    parser_generate.set_defaults(func=generate_clonotypes)

    # simulateReads subparser
    parser_simulate = subparsers.add_parser('simulateReads', help='Simulate sequencing reads')
    parser_simulate.add_argument('--output', type=str, default='reads.fasta', help='Output file for the simulated reads (default: reads.fasta)')
    parser_simulate.add_argument('--clonotypes', '--input', type=str, required=True, help='Input file containing clonotypes')
    parser_simulate.add_argument('--read-length', '--read-len', type=int, default=100, help='Length of each simulated read (default: 100)')
    parser_simulate.add_argument('--read-count', '--num-reads', type=int, default=10000, help='Number of reads to generate (default: 10000)')
    constant_region_group = parser_simulate.add_mutually_exclusive_group()
    constant_region_group.add_argument('--with-constant-region', '--include-constant-region', '--with-c-region', '--include-c-region', dest='with_constant_region', action='store_true', default=True, help='Include reads from the constant region (default)')
    constant_region_group.add_argument('--without-constant-region', '--exclude-constant-region', '--without-c-region', '--exclude-c-region', '--no-c-region', dest='with_constant_region', action='store_false', help='Exclude reads from the constant region')
    parser_simulate.set_defaults(func=simulate_reads)

    # downloadData subparser
    parser_download = subparsers.add_parser('downloadData', parents=[common_parser], help='Download reference data')
    parser_download.add_argument('--output', type=str, help='Output folder (default: data/[species]/[group]/)')
    parser_download.set_defaults(func=download_data_wrapper)

    return parser.parse_args()

def generate_clonotypes(args):
    """
    Handle the generateClonotypes action for multiple groups.

    This function prepares the arguments, checks for the existence of required files for all groups,
    and then calls the simulate_receptor_repertoires function for each specified group.

    Args:
    args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    None
    """
    valid_groups = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]

    if args.group == "all":
        group_list = valid_groups
    else:
        if args.group not in valid_groups:
            print(f"Error: Invalid group '{args.group}'. Valid groups are: {', '.join(valid_groups)} or 'all'.")
            exit(1)
        group_list = [args.group]

    # Check all required files before proceeding
    for group in group_list:
        species_path = f"data/{args.species}/{group}"
        v_file = args.v_file or f"{species_path}/{group}V.fasta"
        j_file = args.j_file or f"{species_path}/{group}J.fasta"
        d_file = args.d_file or f"{species_path}/{group}D.fasta" if receptor_has_d_segment(group) else None
        c_file = args.c_file or f"{species_path}/{group}C.fasta" if args.with_constant_region else None

        try:
            check_file_exists(v_file, 'V')
            check_file_exists(j_file, 'J')
            if d_file:
                check_file_exists(d_file, 'D')
            if c_file:
                check_file_exists(c_file, 'C')
        except FileNotFoundError as e:
            print(f"Error for group {group}: {e}")
            exit(1)  # Exit with error code 1

    # Generate clonotypes for each group
    for group in group_list:
        if args.output is None:
            output_file = f"clonotypes_{group}.fasta"
        else:
            output_file = args.output.replace('.fasta', f'_{group}.fasta')
        
        species_path = f"data/{args.species}/{group}"
        v_file = args.v_file or f"{species_path}/{group}V.fasta"
        j_file = args.j_file or f"{species_path}/{group}J.fasta"
        d_file = args.d_file or f"{species_path}/{group}D.fasta" if receptor_has_d_segment(group) else None
        c_file = args.c_file or f"{species_path}/{group}C.fasta" if args.with_constant_region else None

        simulate_receptor_repertoires(group, output_file, v_file, j_file, d_file, c_file,
                                      num_clonotypes=args.count, shm_rate=args.shm_rate,
                                      apply_shm=args.with_shm, append_constant_region=args.with_constant_region,
                                      preserve_alignment=args.with_alignment_preservation)

def simulate_reads(args):
    """
    Handle the simulateReads action.

    This function calls the generate_reads function to simulate sequencing reads.

    Args:
    args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    None
    """
    try:
        generate_reads(
            input_file=args.clonotypes,
            output_file=args.output,
            read_length=args.read_length,
            read_count=args.read_count,
            with_constant_region=args.with_constant_region
        )
    except FileNotFoundError as e:
        print(f"Error: {str(e)}")
        print("Please check the path to your input file and try again.")
        exit(1)  # Exit with error code 1
    except ValueError as e:
        print(f"Error: {str(e)}")
        print("Please ensure that your input sequences are long enough for the specified read length and other parameters.")
        exit(1)  # Exit with error code 1
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        print("Please check your input and try again. If the problem persists, please report this issue.")
        exit(1)  # Exit with error code 1

def download_data_wrapper(args):
    """
    Handle the downloadData action.

    This function prepares the arguments and calls the download_data function.

    Args:
    args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    None
    """
    if args.output is None:
        args.output = "data"
    
    download_data(args.species, args.group, args.output)

def main():
    """
    Main function to run the AIRR Simulation Tool.

    This function parses command-line arguments and calls the appropriate function
    based on the specified action.

    Returns:
    None
    """
    args = parse_args()
    args.func(args)

if __name__ == "__main__":
    main()