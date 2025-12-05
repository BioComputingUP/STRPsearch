import os, shutil
import sys
import tempfile
import typer
from rich import print as rprint
import mimetypes
from Bio.PDB import PDBParser
pdbparser=PDBParser()
from Bio.PDB import MMCIFParser
cifprarser=MMCIFParser()
import glob
# Add parent directory to sys.path to allow importing from src
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')))

from src import execute_strpsearch as ex
from src import download_structure as ds
from src import config as cfg

# Initialize the Typer app
app = typer.Typer(pretty_exceptions_enable=False, add_completion=False)

# Application version
__version__ = "1.0.0"

@app.command()
def version():
    """
    Show the version of the application and exit.
    """
    rprint(f"[bold green]STRPsearch Lite 2, version: {__version__}[/bold green]\n")
    raise typer.Exit()

def validate_min_height(value: str) -> float:
    """
    Validate that min_height is either "auto" or a numerical value between 0 and 1 (inclusive).
    """
    if value.lower() == "auto":
        return value  # Allow "auto" as a valid value

    try:
        value = float(value)  # Ensure the value is a float
    except ValueError:
        raise typer.BadParameter("'min_height' must be 'auto' or a numerical value between 0 and 1 (inclusive).")

    if value < 0 or value > 1:
        raise typer.BadParameter("'min_height' must be 'auto' or a numerical value between 0 and 1 (inclusive).")
    return value
@app.callback()
def main(
    ctx: typer.Context,
    chain: str = typer.Option("all", help="Specific chain to query from the structures."),
    temp_dir: str = typer.Option(cfg.temp_dir, help="Path to the temporary directory."),
    max_eval: float = typer.Option(cfg.max_eval_p, help="Maximum E-value of the targets to prefilter."),
    min_height: str = typer.Option(cfg.min_height_p, help="Minimum height of TM-score signals to be processed."),
    keep_temp: bool = typer.Option(cfg.keep_temp, help="Whether to keep the temporary directory and files."),
    pymol_pse: bool = typer.Option(cfg.pymol_pse, help="Whether to create and output PyMOL session files."),
    db: str = typer.Option(None, help="Path to the database to use."),
    chainsaw: bool = typer.Option(False, help="Whether to use Chainsaw for domain trimming.")
):
    """
    Main callback to store global options in the context.
    """
    ctx.obj = {
        "chain": chain,
        "temp_dir": temp_dir,
        "max_eval": max_eval,
        "min_height": min_height,
        "keep_temp": keep_temp,
        "pymol_pse": pymol_pse,
        "db": db,
        "chainsaw": chainsaw
    }

@app.command()
def query_file(
    input_file: str = typer.Argument(..., help="Path to the input structure file to query (PDB/mmCIF)."),
    out_dir: str = typer.Argument(..., help="Path to the output directory."),
    chain: str = typer.Option("all", help="Specific chain to query from the structures."),
    temp_dir: str = typer.Option(cfg.temp_dir, help="Path to the temporary directory."),
    max_eval: float = typer.Option(cfg.max_eval_p, help="Maximum E-value of the targets to prefilter." , min=0),
    min_height: str = typer.Option(cfg.min_height_p, help="Minimum height of TM-score signals to be processed.",callback=validate_min_height),
    keep_temp: bool = typer.Option(cfg.keep_temp, help="Whether to keep the temporary directory and files."),
    pymol_pse: bool = typer.Option(cfg.pymol_pse, help="Whether to create and output PyMOL session files."),
    db: str = typer.Option(None, help="Path to the database to use."),
    chainsaw: bool = typer.Option(False, help="Whether to use Chainsaw for domain trimming.")
):
    """
    Query an existing PDB/CIF formatted structure file by providing the file path.
    """
    # Determine database paths
    if db:
        tul_db = os.path.join(db, "tul_foldseek_db", "db")
        rul_db = os.path.join(db, "rul_structure_db")
    else:
        tul_db = os.path.join(cfg.project_root, "data", "databases", "tul_foldseek_db", "db")
        rul_db = os.path.join(cfg.project_root, "data", "databases", "rul_structure_db")
    #verify if query file exists
    if not os.path.exists(input_file):
        rprint(f"[bold red]Input file does not exist: {input_file}[/bold red]")
        sys.exit()

    #get name of input file
    input_filename = os.path.basename(input_file)
    head, tail=os.path.split(input_filename)
    input_file_name=tail.split('.')[0]
    #verify type of query file
    mime_type, encoding = mimetypes.guess_type(input_file)
    if mime_type:
        #if the file if pdb type
        if "pdb" in mime_type or input_file.endswith(".ent.gz") or input_filename.endswith(".ent"):

            pdb_id = ds.extract_structure_and_chains(input_file)[0]#get pdb id from pdb file
            input_name=pdb_id
            out_dir = os.path.join(out_dir, input_name) #construct outdir path as /out_dir/input_name
            print(out_dir)
            #check if output directory exists:
            if os.path.exists(out_dir):
                #if it does then clear it and reuse it
                print(f"alredy exists {out_dir}")
                rprint("[yellow]Warning: Output directory already exists. Reusing it.[/yellow]\n")
                for filen in os.listdir(out_dir):
                    file_path=os.path.join(out_dir,filen)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print('Failed to delete %s. Reason: %s' % (file_path, e))

            else:
                #if not create one
                os.makedirs(out_dir)
            
            #creates temporary directory:
            if not os.path.exists(temp_dir):
                try:
                    os.makedirs(temp_dir)
                    rprint(f"[bold yellow]Created temporary directory: {temp_dir}[/bold yellow]")
                except Exception as e:
                    rprint(f"[bold red]Failed to create temporary directory '{temp_dir}': {e}[/bold red]")
                    sys.exit(1)
            query_dir = os.path.join(out_dir, "query_structures")
            os.makedirs(query_dir, exist_ok=True)
            success= ds.download_pdb_structure(pdb_id=pdb_id, chain=chain, out_dir=query_dir, temp_dir=temp_dir)
            
        elif "cif" in mime_type:
            max_retries = 3
            for attempt in range(1, max_retries + 1):
                try:
                    pdb_id = ds.extract_structure_and_chains_cif(input_file)#get pdb id from pdb file
                    input_name=pdb_id
                    out_dir = os.path.join(out_dir, input_name)
                    if os.path.exists(out_dir):
                        rprint("[yellow]Warning: Output directory already exists. Reusing it.[/yellow]\n")
                        for filen in os.listdir(out_dir):
                            file_path=os.path.join(out_dir,filen)
                            try:
                                if os.path.isfile(file_path) or os.path.islink(file_path):
                                    os.unlink(file_path)
                                elif os.path.isdir(file_path):
                                    shutil.rmtree(file_path)
                            except Exception as e:
                                print('Failed to delete %s. Reason: %s' % (file_path, e))
                    else:
                        os.makedirs(out_dir)
                    if not os.path.exists(temp_dir):
                        try:
                            os.makedirs(temp_dir)
                            rprint(f"[bold yellow]Created temporary directory: {temp_dir}[/bold yellow]")
                        except Exception as e:
                            rprint(f"[bold red]Failed to create temporary directory '{temp_dir}': {e}[/bold red]")
                            sys.exit(1)
                    query_dir = os.path.join(out_dir, "query_structures")
                    os.makedirs(query_dir, exist_ok=True)
                    success , pdb_id, test= ds.extract_chains(input_file=input_file, chain=chain, out_dir=query_dir, temp_dir=temp_dir)
                    
                    break  # Success, exit the retry loop
                except Exception as e:
                    rprint(f"[bold yellow]Attempt {attempt} failed to extract chains: {e}[/bold yellow]")
                    success = False
                    if attempt == max_retries:
                        rprint(f"[bold red]All {max_retries} attempts failed. Exiting.[/bold red]")
                        sys.exit(0)
        else:
            rprint(f"[bold red]Only PDB / mmCIF format is accepted for query files[bold red]\n")
            return False
    else:
        rprint(f"[bold red]The query file format is ambiguous for query {input_file}[bold red]\n")
        return False
    if not success:
        rprint("[bold red]Chain extraction failed.[/bold red]")
        return 
    
    # If a specific chain is requested, verify its file exists
    if chain and chain.lower() != "all":
        chain_file = os.path.join(query_dir, f"{input_file_name}_{chain}.cif")
        if not os.path.isfile(chain_file):
            rprint(f"[bold red]❌ Chain '{chain}' not found in the structure.[/bold red]")
            rprint(f"[bold red]Make sure the chain ID exists in the input file: {input_file}[/bold red]")
            sys.exit(1)

        # Remove any other chain files
        for f in os.listdir(query_dir):
            if f != os.path.basename(chain_file):
                os.remove(os.path.join(query_dir, f))

    # Create a unique subdirectory for this run
    temp_dir = tempfile.mkdtemp(dir=temp_dir)
    rprint(f"[bold green]Using temporary subdirectory: {temp_dir}[/bold green]")

    # Execute the prediction
    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        pymol_pse=pymol_pse,
        max_eval_p=max_eval,
        min_height_p=min_height,
        tul_db=tul_db,
        rul_db=rul_db,
        pdb_id=pdb_id,
        chainsaw=chainsaw)

@app.command()
def download_pdb(
    pdb_id: str = typer.Argument(..., help="PDB ID of the experimental structure to download and query."),
    out_dir: str = typer.Argument(..., help="Path to the output directory."),
    chain: str = typer.Option("all", help="Specific chain to query from the structures."),
    temp_dir: str = typer.Option(cfg.temp_dir, help="Path to the temporary directory."),
    max_eval: float = typer.Option(cfg.max_eval_p, help="Maximum E-value of the targets to prefilter." , min=0),
    min_height: str = typer.Option(cfg.min_height_p, help="Minimum height of TM-score signals to be processed." ,callback=validate_min_height),
    keep_temp: bool = typer.Option(cfg.keep_temp, help="Whether to keep the temporary directory and files."),
    pymol_pse: bool = typer.Option(cfg.pymol_pse, help="Whether to create and output PyMOL session files."),
    db: str = typer.Option(None, help="Path to the database to use."),
    chainsaw: bool = typer.Option(False, help="Whether to use Chainsaw for domain trimming.")
):
    """
    Download and query a structure from the PDB database by providing the PDB ID and the specific chain of interest.
    """
     # Determine database paths
    if db:
        tul_db = os.path.join(db, "tul_foldseek_db", "db")
        rul_db = os.path.join(db, "rul_structure_db")
    else:
        tul_db = os.path.join(cfg.project_root, "data", "databases", "tul_foldseek_db", "db")
        rul_db = os.path.join(cfg.project_root, "data", "databases", "rul_structure_db")
    # Ensure the output directory exists
    if os.path.exists(out_dir):
        rprint("[yellow]Warning: Output directory already exists. Reusing it.[/yellow]\n")
    else:
        os.makedirs(out_dir)

    # Ensure the temporary directory exists
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
            rprint(f"[bold yellow]Created temporary directory: {temp_dir}[/bold yellow]")
        except Exception as e:
            rprint(f"[bold red]Failed to create temporary directory '{temp_dir}': {e}[/bold red]")
            sys.exit(1)





    # Prepare the output directory structure
    out_dir = os.path.join(out_dir, pdb_id)
    query_dir = os.path.join(out_dir, "query_structures")
    os.makedirs(query_dir, exist_ok=True)

    # Download the PDB structure
    success, ch_id = ds.download_pdb_structure(pdb_id=pdb_id, chain=chain, out_dir=query_dir, temp_dir=temp_dir)
    if not success:
        sys.exit()

    # If a specific chain is requested, verify its file exists
    if chain and chain.lower() != "all":
        chain_file = os.path.join(query_dir, f"{pdb_id}_{ch_id}.cif")
        if not os.path.isfile(chain_file):
            rprint(f"[bold red]Chain '{chain}' not found in the {query_dir}.[/bold red]")
            rprint(f"[bold red]Make sure the chain ID exists in the input Protein: {pdb_id}[/bold red]")
            sys.exit(1)

        # Remove any other chain files
        for f in os.listdir(query_dir):
            if f != os.path.basename(chain_file):
                os.remove(os.path.join(query_dir, f))
    
    # Create a unique subdirectory for this run
    temp_dir = tempfile.mkdtemp(dir=temp_dir)
    rprint(f"[bold green]Using temporary subdirectory: {temp_dir}[/bold green]")

    # Execute the prediction
    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        pymol_pse=pymol_pse,
        max_eval_p=max_eval,
        min_height_p=min_height,
        tul_db=tul_db,
        rul_db=rul_db,
        chainsaw=chainsaw,
        pdb_id=pdb_id
    )


@app.command()
def download_model(
    uniprot_id: str = typer.Argument(..., help="UniProt ID of the AlphaFold model."),
    out_dir: str = typer.Argument(..., help="Path to the output directory."),
    version: str = typer.Argument(cfg.af_version, help="AlphaFold model version."),
    temp_dir: str = typer.Option(cfg.temp_dir, help="Path to the temporary directory."),
    max_eval: float = typer.Option(cfg.max_eval_p, help="Maximum E-value of the targets to prefilter." ,min=0),
    min_height: str = typer.Option(cfg.min_height_p, help="Minimum height of TM-score signals to be processed." , callback=validate_min_height),
    keep_temp: bool = typer.Option(cfg.keep_temp, help="Whether to keep the temporary directory and files."),
    pymol_pse: bool = typer.Option(cfg.pymol_pse, help="Whether to create and output PyMOL session files."),
    db: str = typer.Option(None, help="Path to the database to use."),
    chainsaw: bool = typer.Option(False, help="Whether to use Chainsaw for domain trimming.")
):
     # Determine database paths
    if db:
        tul_db = os.path.join(db, "tul_foldseek_db", "db")
        rul_db = os.path.join(db, "rul_structure_db")
    else:
        tul_db = os.path.join(cfg.project_root, "data", "databases", "tul_foldseek_db", "db")
        rul_db = os.path.join(cfg.project_root, "data", "databases", "rul_structure_db")
    """
    Download and query an AlphaFold model by providing the UniProt ID and the AlphaFold version of interest.
    """
    # Ensure the output directory exists
    if os.path.exists(out_dir):
        rprint("[yellow]Warning: Output directory already exists. Reusing it.[/yellow]\n")
    else:
        os.makedirs(out_dir)

    # Ensure the temporary directory exists
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
            rprint(f"[bold yellow]Created temporary directory: {temp_dir}[/bold yellow]")
        except Exception as e:
            rprint(f"[bold red]Failed to create temporary directory '{temp_dir}': {e}[/bold red]")
            sys.exit(1)



    # Prepare the output directory structure
    out_dir = os.path.join(out_dir, uniprot_id)
    query_dir = os.path.join(out_dir, "query_structures")
    os.makedirs(query_dir, exist_ok=True)

    # Download the AlphaFold structure
    success = ds.download_alphafold_structure(uniprot_id, version, query_dir)
    if not success:
        sys.exit()
    
    # Create a unique subdirectory for this run
    temp_dir = tempfile.mkdtemp(dir=temp_dir)
    rprint(f"[bold green]Using temporary subdirectory: {temp_dir}[/bold green]")


    # Execute the prediction
    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        pymol_pse=pymol_pse,
        max_eval_p=max_eval,
        min_height_p=min_height,
        tul_db=tul_db,
        rul_db=rul_db,
        chainsaw=chainsaw,
        pdb_id=uniprot_id
    )


if __name__ == "__main__":
    app()
