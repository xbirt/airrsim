import pytest
import os
import requests
from unittest.mock import patch, mock_open, Mock
from bs4 import BeautifulSoup
from src.download_data import download_segment, download_data

@pytest.fixture
def mock_requests_get():
    with patch('requests.get') as mock_get:
        yield mock_get

def test_download_segment_success(tmp_path):
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = """
    <html>
        <b>Number of results = 1</b>
        <pre>
        >IGHV1-2*01
CAGGTTCAGCTGGTGCAGTCTGGAGCT
        </pre>
    </html>
    """

    output_path = tmp_path / "test_output.fasta"
    
    with patch('requests.get') as mock_requests_get, \
         patch('builtins.open', mock_open()) as mock_file:
        
        mock_requests_get.return_value = mock_response
        
        download_segment("Homo_sapiens", "IGH", "V", str(output_path))

        # Check if the file was opened correctly
        print(f"Debug: mock_file.call_args: {mock_file.call_args}")  # Debug print
        mock_file.assert_called_once_with(str(output_path), 'w')
        
        # Check if the correct content was written
        expected_content = ">IGHV1-2*01\nCAGGTTCAGCTGGTGCAGTCTGGAGCT"
        mock_file().write.assert_called_once_with(expected_content)

    # Verify that requests.get was called with the correct URL
    expected_url = "https://www.imgt.org/genedb/GENElect?query=7.2+IGHV&species=Homo+sapiens"
    mock_requests_get.assert_called_once_with(expected_url)

def test_download_segment_404(mock_requests_get, capsys):
    mock_response = Mock()
    mock_response.status_code = 404
    mock_requests_get.return_value = mock_response

    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to download:" in captured.out

def test_download_segment_no_results(mock_requests_get, capsys):
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = "<html><body>No results found</body></html>"
    mock_requests_get.return_value = mock_response

    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to find result header for:" in captured.out

def test_download_segment_invalid_html(mock_requests_get, capsys):
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = "<html><body>Invalid HTML structure</body></html>"
    mock_requests_get.return_value = mock_response

    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to find result header for:" in captured.out

@patch('src.download_data.download_segment')
@patch('os.makedirs')
def test_download_data_all(mock_makedirs, mock_download_segment):
    download_data("all", "all", "test_output")
    
    assert mock_makedirs.call_count == 14  # 2 species * 7 groups
    assert mock_download_segment.call_count == 48  # 2 species * 7 groups * 4 segments - 2 species * 4 missing D segments

@patch('src.download_data.download_segment')
@patch('os.makedirs')
def test_download_data_specific(mock_makedirs, mock_download_segment):
    download_data("human", "IGH", "test_output")
    
    expected_path = os.path.join("test_output", "human", "IGH")
    mock_makedirs.assert_called_once_with(expected_path, exist_ok=True)
    assert mock_download_segment.call_count == 4  # V, D, J, C

@patch('src.download_data.download_segment')
@patch('os.makedirs')
def test_download_data_no_d_segment(mock_makedirs, mock_download_segment):
    download_data("mouse", "TRA", "test_output")
    
    expected_path = os.path.join("test_output", "mouse", "TRA")
    mock_makedirs.assert_called_once_with(expected_path, exist_ok=True)
    assert mock_download_segment.call_count == 3  # V, J, C (no D)

def test_download_data_invalid_species():
    with pytest.raises(KeyError, match="Invalid species: invalid_species"):
        download_data("invalid_species", "IGH", "test_output")

def test_download_data_invalid_group():
    with pytest.raises(KeyError, match="Invalid group: INVALID"):
        download_data("human", "INVALID", "test_output")

def test_download_data_empty_species():
    # This case is not possible with the current implementation, 
    # but we'll keep the test for future-proofing
    with pytest.raises(KeyError, match="Invalid species: "):
        download_data("", "IGH", "test_output")

def test_download_data_empty_group():
    # This case is not possible with the current implementation, 
    # but we'll keep the test for future-proofing
    with pytest.raises(KeyError, match="Invalid group: "):
        download_data("human", "", "test_output")

@patch('requests.get', side_effect=requests.ConnectionError)
def test_download_segment_connection_error(mock_get, capsys):
    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to download:" in captured.out
    assert "Connection Error" in captured.out

@patch('requests.get')
def test_download_segment_timeout(mock_get, capsys):
    mock_get.side_effect = requests.Timeout
    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to download:" in captured.out

@patch('requests.get')
def test_download_segment_server_error(mock_get, capsys):
    mock_response = Mock()
    mock_response.status_code = 500
    mock_response.text = ""
    mock_get.return_value = mock_response
    download_segment("Homo_sapiens", "IGH", "V", "test_output.fasta")
    captured = capsys.readouterr()
    assert "Failed to download:" in captured.out
    assert "HTTP status: 500" in captured.out