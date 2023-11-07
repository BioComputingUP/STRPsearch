import os
from Bio.PDB import MMCIFParser, PDBIO, PDBList, Select
import requests
import typer
from rich import print


def download_pdb_structure(pdb_id, pdb_chain, out_dir):
    """Downloads the PDB structure of a given PDB ID,
    extracts and saves a given chain at a given output dir"""

    class ChainSelector(Select):

        def __init__(self, target_chain):
            self.target_chain = target_chain

        def accept_chain(self, chain):
            return chain.get_id() == self.target_chain

    # Instantiate essential modules
    cif_parser = MMCIFParser(QUIET=True)
    pdbl = PDBList()

    # Download the PDB structure
    pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format="mmCif")
    temp_pdb_path = os.path.join(out_dir, f"{pdb_id}.cif")
    # Check if the download was successful
    if not os.path.exists(temp_pdb_path):
        print(
            f"[bold red]Failed to download PDB structure {pdb_id}[/bold red]")
        return 1
    else:
        print(f"[bold]PDB structure {pdb_id} was downloaded successfully")

    structure = cif_parser.get_structure(pdb_id, temp_pdb_path)
    os.remove(temp_pdb_path)
    # Create a new structure with only the specified chain
    selector = ChainSelector(pdb_chain)
    io = PDBIO()
    io.set_structure(structure)
    pdb_file_path = os.path.join(out_dir, f"{pdb_id}{pdb_chain}.pdb")
    io.save(pdb_file_path, select=selector)
    print(
        f"[bold]Chain {pdb_chain} of PDB structure {pdb_id} was extracted successfully"
    )


def download_alphafold_structure(uniprot_id, alphafold_version, out_dir):
    """Downloads the AlphaFold structure of a given UniProt ID"""

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        print(f"[bold yellow]Directory {out_dir} already exists[/bold yellow]")

    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{alphafold_version}.pdb"
    response = requests.get(url)
    # Check if the request was successful (status code 200)
    output_path = os.path.join(
        out_dir, f"AF-{uniprot_id}-F1-model_v{alphafold_version}.pdb"
    )
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            file.write(response.content)
        print(f"[bold]{uniprot_id} AlphaFold structure downloaded successfully")
        return response.status_code
    else:
        print(
            f"[bold red]Failed to download {uniprot_id} AlphaFold structure. Status code: {response.status_code}"
        )
        return response.status_code
