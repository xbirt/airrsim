"""
download_data.py

This script handles the downloading of IMGT reference data for adaptive immune receptor repertoire simulation.
It can download data for specific receptor groups and species, or for all available groups and species.

The script downloads V, D (if applicable), J, and C gene segment files from IMGT databases.
"""

import os
import requests
from bs4 import BeautifulSoup
from src.utils import receptor_has_d_segment

def download_segment(species, group, segment, output_path):
    """
    Download the segment data for a specific species, receptor group, and segment type.

    This function handles the downloading of segment data, which is embedded in HTML.

    Args:
    species (str): The species name (either "Homo_sapiens" or "Mus_musculus").
    group (str): The receptor group (e.g., "IGH", "TRA").
    segment (str): The segment type (e.g., "V", "D", "J", "C").
    output_path (str): The local path where the segment file should be saved.

    Returns:
    None
    """
    species_url = "Homo+sapiens" if species == "Homo_sapiens" else "Mus"
    url = f"https://www.imgt.org/genedb/GENElect?query=7.2+{group}{segment}&species={species_url}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # This will raise an HTTPError for bad responses
        
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Find the 'Number of results' header
            result_header = soup.find('b', string=lambda text: text and text.startswith('Number of results ='))
            
            if result_header:
                # Find the next <pre> tag after the result header
                pre_content = result_header.find_next('pre')
                
                if pre_content:
                    # Trim leading and trailing whitespace, then split into lines
                    content_lines = pre_content.text.strip().split('\n')
                    # Remove any remaining blank lines
                    content_lines = [line for line in content_lines if line.strip()]
                    
                    with open(output_path, 'w') as f:
                        f.write('\n'.join(content_lines))
                    print(f"Downloaded: {output_path}")
                else:
                    print(f"Failed to find content after result header for: {url}")
            else:
                print(f"Failed to find result header for: {url}")
        else:
            print(f"Failed to download: {url} (HTTP status: {response.status_code})")
    except requests.ConnectionError:
        print(f"Failed to download: {url} (Connection Error)")
    except requests.HTTPError as e:
        print(f"Failed to download: {url} (HTTP Error: {e.response.status_code})")
    except requests.RequestException as e:
        print(f"Failed to download: {url} ({type(e).__name__})")

def download_data(species, group, output_folder):
    """
    Main function to download data for specified species and groups.

    This function handles the downloading of data for all combinations of specified species and groups.

    Args:
    species (str): The species to download data for. Can be "human", "mouse", or "all".
    group (str): The receptor group to download data for. Can be a specific group or "all".
    output_folder (str): The base folder where all downloaded data should be saved.

    Returns:
    None

    Raises:
    KeyError: If an invalid species or group is provided, or if the resulting species or group list is empty.
    """
    species_map = {
        "human": "Homo_sapiens",
        "mouse": "Mus_musculus"
    }
    valid_groups = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]

    if species == "all":
        species_list = list(species_map.values())
    else:
        species_list = [species_map.get(species)]
        if species_list[0] is None:
            raise KeyError(f"Invalid species: {species}")

    if group == "all":
        group_list = valid_groups
    else:
        if group not in valid_groups:
            raise KeyError(f"Invalid group: {group}")
        group_list = [group]

    for s in species_list:
        species_folder = "human" if s == "Homo_sapiens" else "mouse"
        for g in group_list:
            folder_path = os.path.join(output_folder, species_folder, g)
            os.makedirs(folder_path, exist_ok=True)
            
            # Download V, J, and C files
            for segment in ['V', 'J', 'C']:
                output_path = os.path.join(folder_path, f"{g}{segment}.fasta")
                download_segment(s, g, segment, output_path)
            
            # Download D file if applicable
            if receptor_has_d_segment(g):
                output_path = os.path.join(folder_path, f"{g}D.fasta")
                download_segment(s, g, 'D', output_path)