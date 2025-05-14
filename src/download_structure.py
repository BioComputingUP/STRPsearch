from Bio.PDB import MMCIFParser, PDBParser, Select
from . import general_utils as gu
from rich import print as rprint
import mimetypes
import requests
import os
import gemmi
import gzip
import shutil


class ChainSelector(Select):
    """
    Custom selector class to filter and save specific chains from a structure.
    """

    def __init__(self, target_chain):
        self.target_chain = target_chain

    def accept_chain(self, chain):
        """
        Accept only the target chain for saving.
        """
        return chain.get_id() == self.target_chain and chain.get_parent().id == 0

def extract_structure_and_chains(pdb_file):
    """
    Extracts the PDB ID (structure ID) and chain IDs from a PDB file, including handling .pdb.gz files.

    Args:
        pdb_file (str): Path to the PDB file (can be .pdb or .pdb.gz).

    Returns:
        tuple: A tuple containing the PDB ID (str) and a list of chain IDs (list of str).
    """
    # Handle .pdb.gz files
    decompressed_file = None
    if pdb_file.endswith(".gz"):
        decompressed_file = pdb_file[:-3]  # Remove the .gz extension
        with gzip.open(pdb_file, "rt") as gz_file:  # Open in text mode
            with open(decompressed_file, "w") as out_file:
                out_file.write(gz_file.read())
        pdb_file = decompressed_file  # Update pdb_file to the decompressed file

    # Read the PDB ID from the HEADER line in the file
    pdb_id = None
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HEADER"):
                pdb_id = line[62:66].strip().upper()
                break

    if not pdb_id:
        raise ValueError("PDB ID could not be extracted from the file.")

    # Initialize the PDB parser
    parser = PDBParser(QUIET=True)
    
    # Parse the structure
    structure = parser.get_structure(pdb_id, pdb_file)
    
    # Extract the chain IDs
    chain_ids = [chain.id for chain in structure.get_chains()]
    
    # Clean up decompressed file if it was a .pdb.gz file
    if decompressed_file:
        os.remove(decompressed_file)

    return pdb_id, chain_ids

def extract_chains(input_file, chain, out_dir):
    """
    Extracts specific chains from a PDB/mmCIF file (including .gz compressed files) 
    and saves them as separate PDB or CIF files.

    Args:
        input_file (str): Path to the input structure file (PDB/mmCIF format, optionally .gz compressed).
        chain (str): Chain ID to extract. Use "all" to extract all chains.
        out_dir (str): Directory to save the extracted chain files.

    Returns:
        bool: True if extraction is successful, False otherwise.
    """
    # Handle .gz compressed files
    decompressed_file = None
    if input_file.endswith(".gz"):
        decompressed_file = input_file[:-3]  # Remove the .gz extension
        with gzip.open(input_file, "rb") as gz_file:
            with open(decompressed_file, "wb") as out_file:
                shutil.copyfileobj(gz_file, out_file)
        input_file = decompressed_file  # Update input_file to the decompressed file

    filename = os.path.basename(input_file)[:-4]
    try:
        # Load structure model using gemmi
        doc = gemmi.cif.read_file(input_file)
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)
        model = structure[0]  # First model
        available_chains = {ch.name for ch in model}
    except Exception as e:
        rprint(f"[bold][{gu.time()}][/bold] [bold red]Error reading structure file: {e}[/bold red]")
        return False

    # Determine the file type (PDB or mmCIF) using MIME type
    mime_type, encoding = mimetypes.guess_type(input_file)
    if mime_type:
        if "pdb" in mime_type:
            pdb_parser = PDBParser(QUIET=True)
            structure = pdb_parser.get_structure(filename, input_file)
        elif "cif" in mime_type:
            cif_parser = MMCIFParser(QUIET=True)
            structure = cif_parser.get_structure(filename, input_file)
        else:
            rprint(f"[bold][{gu.time()}][/bold] [bold red]"
                   "Only PDB / mmCIF format is accepted for query files\n")
            return False
    else:
        rprint(f"[bold][{gu.time()}][/bold] [bold red]"
               f"The query file format is ambiguous for query {filename}\n")
        return False

    # Handle chain selection
    if chain == "all":
        chain_list = list(available_chains)
    else:
        chain = chain.upper()
        if chain not in available_chains:
            rprint(f"[bold][{gu.time()}][/bold] [bold red]"
                   f"Chain '{chain}' not found in the structure. Available chains: {', '.join(sorted(available_chains))}\n")
            return False
        chain_list = [chain]
    # Save each selected chain as a separate CIF file
    for ch_id in chain_list:
        new_structure = gemmi.Structure()
        new_model = gemmi.Model(1)
        target_chain = model[ch_id]
        new_model.add_chain(target_chain)
        new_structure.add_model(new_model)

        output_path = os.path.join(out_dir, f"{filename}_{ch_id}.cif")
        new_structure.make_mmcif_document().write_file(output_path)

    # Clean up decompressed file if it was a .gz file
    if decompressed_file:
        os.remove(decompressed_file)

    return True


def download_pdb_structure(pdb_id, out_dir, temp_dir, chain=None):
    """
    Downloads the PDB structure of a given PDB ID and extracts the specified chain(s).

    Args:
        pdb_id (str): PDB ID of the structure to download.
        out_dir (str): Directory to save the extracted chain files.
        temp_dir (str): Temporary directory for downloaded files.
        chain (str, optional): Chain ID to extract. Use "all" to extract all chains.

    Returns:
        bool: True if the download and extraction are successful, False otherwise.
    """
    # Convert PDB ID to lowercase and uppercase for URL compatibility
    pdb_id_lower = pdb_id.lower()
    pdb_id_upper = pdb_id.upper()

    # Create a temporary directory for downloaded structures
    temp_structure_dir = os.path.join(temp_dir, "downloaded_structures")
    os.makedirs(temp_structure_dir, exist_ok=True)

    # Construct the URL for the mmCIF file
    url = f"https://files.rcsb.org/download/{pdb_id_upper}.cif"
    temp_structure_path = os.path.join(temp_structure_dir, f"{pdb_id_lower}.cif")

    # Download the structure using requests
    response = requests.get(url)
    if response.status_code != 200:
        rprint(f"[bold][{gu.time()}][/bold] [bold red]"
               f"Failed to download PDB structure {pdb_id_lower}")
        return False
    else:
        with open(temp_structure_path, "wb") as f:
            f.write(response.content)
        rprint(f"[bold][{gu.time()}][/bold] [bold]"
               f"PDB structure {pdb_id_lower} was downloaded successfully")

    # Extract the desired chain(s) using the extract_chains function
    success = extract_chains(
        input_file=temp_structure_path,
        chain=chain,
        out_dir=out_dir
    )

    return success


def download_alphafold_structure(uniprot_id, alphafold_version, out_dir):
    """
    Downloads the AlphaFold structure of a given UniProt ID.

    Args:
        uniprot_id (str): UniProt ID of the structure to download.
        alphafold_version (str): Version of the AlphaFold model.
        out_dir (str): Directory to save the downloaded structure.

    Returns:
        bool: True if the download is successful, False otherwise.
    """
    # Convert UniProt ID to uppercase for URL compatibility
    uniprot_id = uniprot_id.upper()

    # Construct the URL for the AlphaFold model
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{alphafold_version}.cif"
    rprint(f"[bold][{gu.time()}][/bold] [bold]"
           f"Downloading AlphaFold-predicted model of {uniprot_id}\n")

    # Download the structure using requests
    response = requests.get(url)
    output_path = os.path.join(out_dir, f"AF-{uniprot_id}-F1-model_v{alphafold_version}_A.cif")
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            file.write(response.content)
        rprint(f"[bold][{gu.time()}][/bold] [bold]"
               f"AlphaFold-predicted model of {uniprot_id} was downloaded successfully")
        return True
    else:
        rprint(f"[bold][{gu.time()}][/bold] [bold red]"
               f"Failed to download AlphaFold-predicted model of {uniprot_id}. Status code: {response.status_code}")
        return False
