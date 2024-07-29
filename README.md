# AIRRSIM: Adaptive Immune Receptor Repertoire Simulator

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Installation](#installation)
4. [Usage](#usage)
   - [Generating Clonotypes](#generating-clonotypes)
   - [Simulating Sequencing Reads](#simulating-sequencing-reads)
   - [Downloading Reference Data](#downloading-reference-data)
5. [Detailed Functionality](#detailed-functionality)
6. [File Descriptions](#file-descriptions)
7. [Dependencies](#dependencies)
8. [Contributing](#contributing)
9. [License](#license)

## Introduction

AIRRSIM is a powerful tool designed for simulating Adaptive Immune Receptor Repertoire (AIRR) clonotypes and sequence read data. It provides a comprehensive suite of functionalities for generating clonotypes, simulating sequencing reads, and downloading necessary reference data. This tool is particularly useful for researchers and bioinformaticians working in the field of immunogenomics and adaptive immunity.

## Features

- **Clonotype Generation**: Simulate diverse repertoires of adaptive immune receptors, including immunoglobulins and T-cell receptors.
- **Read Simulation**: Generate synthetic sequencing reads from the simulated clonotypes.
- **Reference Data Download**: Easily obtain required reference data from IMGT databases.
- **Customizable Parameters**: Fine-tune your simulations with a wide range of adjustable parameters.
- **Support for Multiple Species**: Generate data for human and mouse repertoires.
- **Somatic Hypermutation Simulation**: Model the process of somatic hypermutation in B-cell receptors.

## Installation

To install AIRRSIM, follow these steps:

1. Clone the repository:
   
    git clone https://github.com/xbirt/airrsim.git
    cd airrsim

2. Install the required dependencies:
   
    pip install -r requirements.txt

3. Install the package:
   
    pip install .

## Usage

AIRRSIM provides a command-line interface for its main functionalities. Here are the primary commands:

### Generating Clonotypes

Generate clonotypes for a specific receptor group:

    airrsim generateClonotypes <group> [options]

Options:
- `<group>`: Receptor group (TRA, TRB, TRG, TRD, IGH, IGK, IGL, all)
- `--species`: Species (human, mouse, all) [default: human]
- `--output`: Output file [default: clonotypes[group].fasta]
- `--count`: Number of clonotypes to generate [default: 1000]
- `--shm_rate`: Somatic hypermutation rate [default: 0.001]
- `--with-shm`: Apply somatic hypermutations [default: True]
- `--no-shm`: Do not apply somatic hypermutations
- `--include_constant`: Include constant region in the output sequences [default: True]
- `--v_file`, `--d_file`, `--j_file`, `--c_file`: Paths to gene segment files

Example:

    airrsim generateClonotypes IGH --species human --count 10000 --shm_rate 0.002 --output my_igh_clonotypes.fasta

### Simulating Sequencing Reads

Generate simulated sequencing reads from clonotypes:

    airrsim simulateReads [options]

Options:
- `--output`: Output file for the simulated reads [default: reads.fasta]
- `--clonotypes`: Input file containing clonotypes [required]
- `--read-length`: Length of each simulated read [default: 100]
- `--read-count`: Number of reads to generate [default: 10000]
- `--no-c-region`: Avoid generating reads that are located mostly within the constant region

Example:

    airrsim simulateReads --clonotypes my_igh_clonotypes.fasta --read-length 150 --read-count 100000 --output my_simulated_reads.fasta

### Downloading Reference Data

Download reference data for gene segments:

    airrsim downloadData <group> [options]

Options:
- `<group>`: Receptor group (TRA, TRB, TRG, TRD, IGH, IGK, IGL, all)
- `--species`: Species (human, mouse, all) [default: human]
- `--output`: Output folder [default: data/[species]/[group]/]

Example:

    airrsim downloadData IGH --species human --output my_reference_data

## Detailed Functionality

### Clonotype Generation
- Simulates V(D)J recombination process for adaptive immune receptors.
- Incorporates realistic biological features such as P and N nucleotide additions and exonuclease trimming.
- Models somatic hypermutation (SHM) for B-cell receptors with customizable mutation rates.
- Generates diverse repertoires based on input gene segment files.

### Read Simulation
- Creates synthetic sequencing reads from generated clonotypes.
- Allows specification of read length and count.
- Option to avoid generating reads primarily from constant regions.
- Simulates realistic sequencing scenarios for repertoire analysis.

### Reference Data Download
- Fetches V, D, J, and C gene segment files from IMGT databases.
- Supports data download for human and mouse species.
- Organizes downloaded data into a structured directory format.

## File Descriptions

- `src/main.py`: Main entry point for the AIRRSIM tool.
- `src/simulation.py`: Contains functions for simulating receptor repertoires.
- `src/read_simulation.py`: Handles the simulation of sequencing reads.
- `src/download_data.py`: Manages the downloading of reference data.
- `src/recombination.py`: Implements the V(D)J recombination process.
- `src/utils.py`: Utility functions for file handling and sequence processing.

## Dependencies

AIRRSIM requires the following Python packages:

- biopython>=1.79
- numpy>=1.21.0
- requests>=2.26.0
- beautifulsoup4>=4.9.3

These dependencies are listed in the `requirements.txt` file and will be installed automatically when following the installation instructions.

## Contributing

Contributions to AIRRSIM are welcome! Please feel free to submit pull requests, create issues or suggest improvements.

## License

AIRRSIM is released under the MIT License. See the [LICENSE](LICENSE) file for details.

---

For more information or support, please open an issue on the GitHub repository.