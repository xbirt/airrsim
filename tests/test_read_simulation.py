import pytest
import os
from src.read_simulation import parse_sequence_id, generate_reads
from src.utils import read_fasta_with_error_handling

@pytest.fixture
def sample_fasta(tmp_path):
    fasta_content = (
        ">TRB_clonotype_1_TRBV7-9*05_TRBD2*02_TRBJ1-6*01_TRBC1*03_V288-13+2+12_D16-0-4_J1+1+53-0_C108\n"
        "GATACTGGAGTCTCCCAGAACCCCAGACACAAGATCACAAAGAGGGGACAGAATGTAACTTTCAGGTGTGATCCAATTTCTGAACACAACCGCCTTTATTGGTACCGACAGACCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCCAGAATGAAGCTCAACTAGAAAAATCAAGGCTGCTCAGTGATCGGTTCTCTGCAGAGAGGCCTAAGGGATCTCTCTCCACCTTGGAGATCCAGCGCACAGAGCAGGGGGACTCGGCCATGTATCTCTGTCCTAACCCTGGGCGGGACTAGCGGGCCCGCTCCTATAATTCACCCCTCCACTTTGGGAATGGGACCAGGCTCACTGTGACAGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATG\n"
        ">TRB_clonotype_2_TRBV4-1*01_TRBD1*01_TRBJ1-6*01_TRBC2*01_V287-7+2+1_D12-5-5_J11+1+53-13_C18\n"
        "GACACTGAAGTTACCCAGACACCAAAACACCTGGTCATGGGAATGACAAATAAGAAGTCTTTGAAATGTGAACAACATATGGGGCACAGGGCTATGTATTGGTACAAGCAGAAAGCTAAGAAGCCACCGGAGCTCATGTTTGTCTACAGCTATGAGAAACTCTCTATAAATGAAAGTGTGCCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCTCTTAAACCTTCACCTACACGCCCTGCAGCCAGAAGACTCAGCCCTGTATCTCTGCGCCAGCATCCCAGTCGTCAGACATTGCCCCTCCACTTTGGGAATGGGACCAGGCTCACTGTGACAGGACTGTGGCTTCACCTCC\n"
        ">TRB_clonotype_3_TRBV6-4*01_TRBD2*02_TRBJ2-1*01_TRBC2*02_V287-5+0+10_D16-4-3_J3+2+50-7_C24\n"
        "ATTGCTGGGATCACCCAGGCACCAACATCTCAGATCCTGGCAGCAGGACGGCGCATGACACTGAGATGTACCCAGGATATGAGACATAATGCCATGTACTGGTATAGACAAGATCTAGGACTGGGGCTAAGGCTCATCCATTATTCAAATACTGCAGGTACCACTGGCAAAGGAGAAGTCCCTGATGGTTATAGTGTCTCCAGAGCAAACACAGATGATTTCCCCCTCACGTTGGCGTCTGCTGTACCCTCTCAGACATCTGTGTACTTCTGTGCCAGCAGTTCTCTACTCTTCTAGCGGGAGAGTTTAATGAGCAGTTCTTCGGGCCAGGGACACGGCTCACCGTGCTAGGTCAAGAGAAAGGATTCCAGAGGC\n"
    )
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)

def test_parse_sequence_id():
    assert parse_sequence_id("TRB_clonotype_1_TRBV7-9*05_TRBD2*02_TRBJ1-6*01_TRBC1*03_V288-13+2+12_D16-0-4_J1+1+53-0_C108") == (55, 108)
    assert parse_sequence_id("TRB_clonotype_2_TRBV4-1*01_TRBD1*01_TRBJ1-6*01_TRBC2*01_V287-7+2+1_D12-5-5_J11+1+53-13_C18") == (52, 18)
    assert parse_sequence_id("TRB_clonotype_3_TRBV6-4*01_TRBD2*02_TRBJ2-1*01_TRBC2*02_V287-5+0+10_D16-4-3_J3+2+50-7_C24") == (48, 24)
    assert parse_sequence_id("Invalid_ID") == (None, None)

def test_generate_reads(sample_fasta, tmp_path):
    output_file = str(tmp_path / "output.fasta")
    generate_reads(sample_fasta, output_file, read_length=50, read_count=100)
    
    assert os.path.exists(output_file)
    reads = read_fasta_with_error_handling(output_file, "test")
    
    assert len(reads) == 100
    for read_id, read_seq in reads:
        assert len(read_seq) == 50
        assert read_id.startswith(("TRB_c1_R", "TRB_c2_R", "TRB_c3_R"))

def test_generate_reads_no_c_region(sample_fasta, tmp_path):
    output_file = str(tmp_path / "output_no_c.fasta")
    generate_reads(sample_fasta, output_file, read_length=50, read_count=100, no_c_region=True)
    
    assert os.path.exists(output_file)
    reads = read_fasta_with_error_handling(output_file, "test")
    
    # Read the original clonotypes
    original_clonotypes = read_fasta_with_error_handling(sample_fasta, "Input clonotype")
    clonotype_dict = {f"c{i+1}": seq for i, (_, seq) in enumerate(original_clonotypes)}
    
    assert len(reads) == 100
    for read_id, read_seq in reads:
        assert len(read_seq) == 50
        assert read_id.startswith(("TRB_c1_R", "TRB_c2_R", "TRB_c3_R"))
        
        clonotype_id = read_id.split('_')[1]  # e.g., 'c1', 'c2', 'c3'
        start_pos = int(read_id.split('_R')[1].split('_')[0])
        j_length, c_length = parse_sequence_id(read_id)
        
        original_seq = clonotype_dict[clonotype_id]
        assert start_pos <= len(original_seq) - c_length - j_length - 12 - len(read_seq)
        
        # Additional check to ensure the read actually matches the original sequence
        assert read_seq == original_seq[start_pos:start_pos+len(read_seq)]

def test_generate_reads_file_not_found():
    with pytest.raises(FileNotFoundError):
        generate_reads("nonexistent.fasta", "output.fasta", read_length=50, read_count=100)

def test_generate_reads_invalid_input(tmp_path):
    # Test case 1: Sequence shorter than read length
    short_seq_fasta = tmp_path / "short_seq.fasta"
    short_seq_fasta.write_text(">TRB_clonotype_1_V1_J1_C1\nACGT\n")
    output_file = str(tmp_path / "output_short.fasta")
    
    with pytest.raises(ValueError, match="Sequence length .* is smaller than read length"):
        generate_reads(str(short_seq_fasta), output_file, read_length=50, read_count=100)
    
    # The output file might be created but should be empty
    assert os.path.getsize(output_file) == 0  # Assuming an empty or nearly empty file

    # Test case 2: Sequence too short for no_c_region option
    short_seq_no_c_fasta = tmp_path / "short_seq_no_c.fasta"
    short_seq_no_c_fasta.write_text(">TRB_clonotype_1_V1_J10_C10\n" + "A" * 60 + "\n")
    output_file_no_c = str(tmp_path / "output_short_no_c.fasta")
    
    with pytest.raises(ValueError, match="Sequence is too short to generate reads with the given parameters"):
        generate_reads(str(short_seq_no_c_fasta), output_file_no_c, read_length=50, read_count=100, no_c_region=True)
    
    # The output file might be created but should be empty
    assert os.path.getsize(output_file_no_c) == 0

    # Test case 3: Invalid FASTA format
    invalid_fasta = tmp_path / "invalid.fasta"
    invalid_fasta.write_text("Invalid FASTA content")
    output_file_invalid = str(tmp_path / "output_invalid.fasta")
    
    generate_reads(str(invalid_fasta), output_file_invalid, read_length=50, read_count=100)
    
    # In this case, the output file should not be created
    assert not os.path.exists(output_file_invalid)

    # Test case 4: Non-existent input file
    non_existent_file = str(tmp_path / "non_existent.fasta")
    output_file_non_existent = str(tmp_path / "output_non_existent.fasta")
    
    with pytest.raises(FileNotFoundError, match="Input file not found"):
        generate_reads(non_existent_file, output_file_non_existent, read_length=50, read_count=100)
    
    # In this case, the output file should not be created
    assert not os.path.exists(output_file_non_existent)

@pytest.mark.parametrize("read_length,read_count", [
    (30, 50),
    (100, 200),
    (200, 10),
])
def test_generate_reads_different_params(sample_fasta, tmp_path, read_length, read_count):
    output_file = str(tmp_path / f"output_{read_length}_{read_count}.fasta")
    generate_reads(sample_fasta, output_file, read_length=read_length, read_count=read_count)
    
    assert os.path.exists(output_file)
    reads = read_fasta_with_error_handling(output_file, "test")
    
    assert len(reads) == read_count
    for read_id, read_seq in reads:
        assert len(read_seq) == read_length

def test_generate_reads_incompatible_id_format(tmp_path):
    fasta_content = (
        ">Incompatible_ID_1\n"
        "CAAGGAGGTGGAGCAGGATCCTCAGGCCTTGAGACAGACCCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGACTTGGTCAAGAGAAATGACATCCCTAAACATCTATGGCCATGAGAAACTCTTTGGCATCGGCCTCTTCCTCCTCCTCCTCAATGCTGCTTGCACATACTTCTGTGCCACCAGTATAGACATCCT\n"
    )
    fasta_file = tmp_path / "incompatible.fasta"
    fasta_file.write_text(fasta_content)
    output_file = str(tmp_path / "output.fasta")
    
    with pytest.warns(UserWarning, match="Cannot avoid constant region due to incompatible sequence ID format"):
        generate_reads(str(fasta_file), output_file, read_length=50, read_count=100, no_c_region=True)
    
    assert os.path.exists(output_file)
    reads = read_fasta_with_error_handling(output_file, "test")
    assert len(reads) == 100