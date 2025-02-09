from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from . import general_utils as gu
from rich import print as rprint
import mimetypes
import requests
import os


class ChainSelector(Select):

    def __init__(self, target_chain):
        self.target_chain = target_chain

    def accept_chain(self, chain):
        if chain.get_id() == self.target_chain and chain.get_parent().id == 0:
            return True


def extract_chains(input_file, chain, out_dir):
    filename = os.path.basename(input_file)[:-4]

    mime_type, encoding = mimetypes.guess_type(input_file)
    if mime_type:
        if "pdb" in mime_type:
            pdb_parser = PDBParser(QUIET=True)
            # Parse the structure from the PDB file
            structure = pdb_parser.get_structure(filename, input_file)
        elif "cif" in mime_type:
            # Instantiate essential modules
            cif_parser = MMCIFParser(QUIET=True)
            # Parse the structure from the PDB file
            structure = cif_parser.get_structure(filename, input_file)
        else:
            rprint(f"[bold][{gu.time()}][bold] [bold red]"
                   "Only PDB / mmCIF format is accepted for query files\n")
            # logging.error("Only PDB / mmCIF format is accepted for query files")
            return False
    else:
        rprint(f"[bold][{gu.time()}][bold] [bold red]"
               f"The query file format is ambiguous for query {filename}\n")
        # logging.error(f"The query file format is ambiguous for query {query_name}")
        return False

    structure_model = structure[0]
    chain_list = []
    if chain == "all":
        for chain in structure_model:
            chain_list.append(chain.id)
    else:
        chain = chain.upper()
        chain_list.append(chain)

    for chain in chain_list:
        # Create a new structure with only the specified chain
        selector = ChainSelector(chain)
        io = PDBIO()
        io.set_structure(structure)
        output_path = os.path.join(out_dir, f"{filename}_{chain}.pdb")
        io.save(output_path, select=selector)
    return True


def download_pdb_structure(pdb_id, out_dir, temp_dir, chain=None):
    """
    Downloads the PDB structure of a given PDB ID.
    Extracts and saves a given chain at the specified output directory.
    Returns the success status.
    """
    # Store pdb_id in both lowercase and uppercase forms.
    # RCSB URLs require the PDB code in uppercase.
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

    # Extract the desired chain(s) using your existing function
    success = extract_chains(
        input_file=temp_structure_path,
        chain=chain,
        out_dir=out_dir
    )

    return success


def download_alphafold_structure(uniprot_id, alphafold_version, out_dir):
    """
    Downloads the AlphaFold structure of a given UniProt ID
    Return the success status
    """

    uniprot_id = uniprot_id.upper()

    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{alphafold_version}.pdb"
    rprint(f"[bold][{gu.time()}][/bold] [bold]"
           f"Downloading AlphaFold-predicted model of {uniprot_id}\n")
    response = requests.get(url)
    # Check if the request was successful (status code 200)
    output_path = os.path.join(out_dir, f"AF-{uniprot_id}-F1-model_v{alphafold_version}_A.pdb")
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
