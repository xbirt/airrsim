import pytest
import os
from src.utils import receptor_has_d_segment, read_fasta_with_error_handling, read_fasta, save_fasta

# Test data
TEST_FASTA_CONTENT = """>seq1
ATCG
>seq2
GCTA
"""

@pytest.fixture
def temp_fasta_file(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(TEST_FASTA_CONTENT)
    return str(fasta_file)

def test_receptor_has_d_segment():
    assert receptor_has_d_segment('IGH') == True
    assert receptor_has_d_segment('TRB') == True
    assert receptor_has_d_segment('TRD') == True
    assert receptor_has_d_segment('TRA') == False
    assert receptor_has_d_segment('IGK') == False

def test_read_fasta_with_error_handling(temp_fasta_file):
    result = read_fasta_with_error_handling(temp_fasta_file, "V")
    assert len(result) == 2
    assert result[0] == ('seq1', 'ATCG')
    assert result[1] == ('seq2', 'GCTA')

def test_read_fasta_with_error_handling_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_fasta_with_error_handling("nonexistent_file.fasta", "V")

def test_read_fasta_with_error_handling_empty_file(tmp_path):
    empty_file = tmp_path / "empty.fasta"
    empty_file.write_text("")
    with pytest.raises(ValueError):
        read_fasta_with_error_handling(str(empty_file), "V")

def test_read_fasta(temp_fasta_file):
    result = list(read_fasta(temp_fasta_file))
    assert len(result) == 2
    assert result[0] == ('seq1', 'ATCG')
    assert result[1] == ('seq2', 'GCTA')

def test_read_fasta_nonexistent_file(capfd):
    list(read_fasta("nonexistent_file.fasta"))
    captured = capfd.readouterr()
    assert "Warning: File nonexistent_file.fasta does not exist." in captured.out

def test_save_fasta(tmp_path):
    output_file = tmp_path / "output.fasta"
    sequences = [('seq1', 'ATCG' * 20), ('seq2', 'GCTA' * 25)]
    save_fasta(sequences, str(output_file), line_length=40)
    
    assert output_file.exists()
    content = output_file.read_text()
    assert ">seq1" in content
    assert ">seq2" in content
    assert "ATCGATCGATCGATCGATCGATCGATCGATCGATCG" in content
    assert "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" in content

def test_save_fasta_line_length(tmp_path):
    output_file = tmp_path / "output_line_length.fasta"
    sequences = [('seq1', 'A' * 100)]
    save_fasta(sequences, str(output_file), line_length=50)
    
    content = output_file.read_text().split('\n')
    assert len(content[1]) == 50
    assert len(content[2]) == 50