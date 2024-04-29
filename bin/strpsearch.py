import os
import sys

# Add parent directory to sys.path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src import execute_strpsearch as ex
from src import download_structure as ds
from typing_extensions import Annotated
from rich import print as rprint
from src import config as cfg
import tempfile
import typer


# Initialize the app
# Disable pretty_exceptions
app = typer.Typer(pretty_exceptions_enable=False, add_completion=False)

__version__ = "0.1.0"


@app.command()
def version():
    """
    Show the version and exit.
    """
    rprint(f"[bold green]RepeatsDB Lite 2, version: {__version__}[/bold green]\n")
    raise typer.Exit()


@app.command()
def query_file(
        input_file: Annotated[
            str, typer.Argument(
                help="Path to the input structure file to query (PDB/mmCIF)"
            )
        ],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to the output directory"
            )
        ],
        chain: Annotated[
            str, typer.Option(
                help="Specific chain to query from the structure"
            )
        ] = cfg.chain,
        temp_dir: Annotated[
            str, typer.Option(
                help="Path to the temporary directory"
            )
        ] = cfg.temp_dir,
        keep_temp: Annotated[
            bool, typer.Option(
                help="Whether to keep the temporary directory and files"
            )
        ] = cfg.keep_temp,
        max_eval: Annotated[
            float, typer.Option(
                help="Maximum E-value of the hits to prefilter"
            )
        ] = cfg.max_eval_p,
        min_height: Annotated[
            str, typer.Option(
                help="Minimum height of TM-score signals to be processed"
            )
        ] = cfg.min_height_p,
):
    """
    Query an existing PDB/CIF formatted structure file by providing the file path.
    """

    if os.path.exists(out_dir):
        rprint("[bold red]Output directory already exists[/bold red]\n")
        sys.exit()
    else:
        os.makedirs(out_dir)

    if not os.path.exists(temp_dir):
        rprint("[bold red]Temporary directory does not exist[/bold red]\n")
        sys.exit()
    else:
        temp_dir = tempfile.mkdtemp(dir=temp_dir)

    if not os.path.exists(input_file):
        rprint(f"[bold red]Input file does not exist")
        sys.exit()

    if min_height != "auto":
        try:
            min_height = float(min_height)
        except ValueError:
            rprint("[bold red]'min_height' parameter must be a numerical value between 0 or 1[/bold red]\n")
            sys.exit()
        if min_height < 0 or min_height > 1:
            rprint("[bold red]'min_height' parameter cannot be smaller than 0 or bigger than 1[/bold red]\n")
            sys.exit()

    query_dir = os.path.join(out_dir, "query_structures")
    os.makedirs(query_dir)

    success = ds.extract_chains(input_file=input_file, chain=chain, out_dir=query_dir)
    if not success:
        sys.exit()

    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        max_eval_p=max_eval,
        min_height_p=min_height
    )


@app.command()
def download_pdb(
        pdb_id: Annotated[
            str, typer.Argument(
                help="PDB ID of the experimental structure to download and query"
            )
        ],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to the output directory"
            )
        ],
        chain: Annotated[
            str, typer.Option(
                help="Specific chain to query from the structure"
            )
        ] = cfg.chain,
        temp_dir: Annotated[
            str, typer.Option(
                help="Path to the temporary directory"
            )
        ] = cfg.temp_dir,
        keep_temp: Annotated[
            bool, typer.Option(
                help="Whether to keep the temporary directory and files"
            )
        ] = cfg.keep_temp,
        max_eval: Annotated[
            float, typer.Option(
                help="Maximum E-value of the hits to prefilter"
            )
        ] = cfg.max_eval_p,
        min_height: Annotated[
            str, typer.Option(
                help="Minimum height of TM-score signals to be processed"
            )
        ] = cfg.min_height_p
):
    """
    Download and query a structure from PDB by providing the PDB ID and the specific Chain of interest.
    """

    if os.path.exists(out_dir):
        rprint("[bold red]Output directory already exists[/bold red]\n")
        sys.exit()
    else:
        os.makedirs(out_dir)

    if not os.path.exists(temp_dir):
        rprint("[bold red]Temporary directory does not exist[/bold red]\n")
        sys.exit()
    else:
        temp_dir = tempfile.mkdtemp(dir=temp_dir)

    if min_height != "auto":
        try:
            min_height = float(min_height)
        except ValueError:
            rprint("[bold red]'min_height' parameter must be a numerical value between 0 or 1[/bold red]\n")
            sys.exit()
        if min_height < 0 or min_height > 1:
            rprint("[bold red]'min_height' parameter cannot be smaller than 0 or bigger than 1[/bold red]\n")
            sys.exit()

    query_dir = os.path.join(out_dir, "query_structures")
    os.makedirs(query_dir)

    success = ds.download_pdb_structure(pdb_id=pdb_id, chain=chain, out_dir=query_dir, temp_dir=temp_dir)
    if not success:
        sys.exit()

    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        max_eval_p=max_eval,
        min_height_p=min_height
    )


@app.command()
def download_model(
        uniprot_id: Annotated[
            str, typer.Argument(
                help="UniProt ID of the AlphaFold-predicted model to download and query"
            )
        ],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to the output directory"
            )
        ],
        af_version: Annotated[
            str, typer.Option(
                help="Version of AlphaFold to download predicted models from"
            )
        ] = cfg.af_version,
        temp_dir: Annotated[
            str, typer.Option(
                help="Path to the temporary directory"
            )
        ] = cfg.temp_dir,
        keep_temp: Annotated[
            bool, typer.Option(
                help="Whether to keep the temporary directory and files"
            )
        ] = cfg.keep_temp,
        max_eval: Annotated[
            float, typer.Option(
                help="Maximum E-value of the hits to prefilter"
            )
        ] = cfg.max_eval_p,
        min_height: Annotated[
            str, typer.Option(
                help="Minimum height of TM-score signals to be processed"
            )
        ] = cfg.min_height_p
):
    """
    Download and query an AlphaFold model by providing the UniProt ID and the AlphaFold version of interest.
    """

    if os.path.exists(out_dir):
        rprint("[bold red]Output directory already exists[/bold red]\n")
        sys.exit()
    else:
        os.makedirs(out_dir)

    if not os.path.exists(temp_dir):
        rprint("[bold red]Temporary directory does not exist[/bold red]\n")
        sys.exit()
    else:
        temp_dir = tempfile.mkdtemp(dir=temp_dir)

    if min_height != "auto":
        try:
            min_height = float(min_height)
        except ValueError:
            rprint("[bold red]'min_height' parameter must be a numerical value between 0 or 1[/bold red]\n")
            sys.exit()
        if min_height < 0 or min_height > 1:
            rprint("[bold red]'min_height' parameter cannot be smaller than 0 or bigger than 1[/bold red]\n")
            sys.exit()

    query_dir = os.path.join(out_dir, "query_structures")
    os.makedirs(query_dir)

    success = ds.download_alphafold_structure(uniprot_id, af_version, query_dir)
    if not success:
        sys.exit()

    ex.execute_predstrp(
        structure_dir=query_dir,
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
        max_eval_p=max_eval,
        min_height_p=min_height
    )


if __name__ == "__main__":
    app()
