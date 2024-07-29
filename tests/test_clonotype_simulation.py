import pytest
import os
from unittest.mock import patch, mock_open, MagicMock
from src.clonotype_simulation import (
    extract_segment_id,
    load_constant_regions,
    get_constant_region,
    generate_clonotypes,
    simulate_receptor_repertoires
)

# Test data
MOCK_FASTA_DATA = (
    ">IGHV1-18*01|IGHV1-18*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+0=296| | |\n"
    "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC\n"
    "TCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCC\n"
    "CCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTAT\n"
    "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTAC\n"
    "ATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGA\n"
)

@pytest.fixture
def mock_fasta_file(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(MOCK_FASTA_DATA)
    return str(fasta_file)

def test_extract_segment_id():
    header = ">IGHV1-18*01|IGHV1-18*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+0=296| | |"
    assert extract_segment_id(header) == "IGHV1-18*01"

    header_no_match = ">Invalid header format"
    assert extract_segment_id(header_no_match) == "unknown"

def test_load_constant_regions(mock_fasta_file):
    constant_regions = load_constant_regions(mock_fasta_file)
    assert len(constant_regions) == 1
    assert constant_regions[0][0].startswith("IGHV1-18*01")
    assert constant_regions[0][1] == "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGA"

def test_load_constant_regions_empty_file(tmp_path):
    empty_file = tmp_path / "empty.fasta"
    empty_file.write_text("")
    with pytest.raises(ValueError, match="No constant region sequences found in file"):
        load_constant_regions(str(empty_file))

def test_load_constant_regions_invalid_file():
    with pytest.raises(ValueError, match="Error reading the constant region file"):
        load_constant_regions("nonexistent_file.fasta")

def test_get_constant_region():
    mock_constant_regions = [
        (">X00000|IGHM|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEM"),
        (">X00000|IGHG|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEG"),
        (">X00000|IGHA|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEA"),
        (">X00000|IGHD|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCED"),
        (">X00000|IGHE|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEE"),
        (">X00000|IGLC|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEIGL"),
        (">X00000|IGKC|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCEIGK"),
        (">X00000|TRBC|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCETRBC"),
        (">X00000|TRAC|Homo sapiens|F|CH1|1..294|294 nt|1| | | | |294+0=294| | |", "SEQUENCETRAC")
    ]

    # Test IGH constant region selection
    with patch('random.choices', return_value=['IGHM']):
        const_seq, const_id = get_constant_region("IGH", mock_constant_regions)
        assert const_seq == "SEQUENCEM"
        assert const_id == "IGHM"

    # Test IGL constant region selection
    const_seq, const_id = get_constant_region("IGL", mock_constant_regions)
    assert const_seq == "SEQUENCEIGL"
    assert const_id == "IGLC"

    # Test IGK constant region selection
    const_seq, const_id = get_constant_region("IGK", mock_constant_regions)
    assert const_seq == "SEQUENCEIGK"
    assert const_id == "IGKC"

    # Test TRB constant region selection
    const_seq, const_id = get_constant_region("TRB", mock_constant_regions)
    assert const_seq == "SEQUENCETRBC"
    assert const_id == "TRBC"

    # Test TRA constant region selection
    const_seq, const_id = get_constant_region("TRA", mock_constant_regions)
    assert const_seq == "SEQUENCETRAC"
    assert const_id == "TRAC"

    # Test unknown receptor type
    with pytest.raises(ValueError, match="Unknown receptor type: UNKNOWN"):
        get_constant_region("UNKNOWN", mock_constant_regions)

    # Test no matching constant region
    with pytest.raises(ValueError, match="No matching constant region found for TRD"):
        get_constant_region("TRD", mock_constant_regions)

    # Test empty constant regions list
    with pytest.raises(ValueError, match="No matching constant region found for IGH"):
        get_constant_region("IGH", [])