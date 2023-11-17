import os
import sys

# Add parent directory to sys.path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src import config as cfg
import typer
from typing_extensions import Annotated
from src import execute_repeatsalgorithm as ex
from rich import print as rprint
from src import download_structure as ds

# Initialize the app
# Disable pretty_exceptions
app = typer.Typer(pretty_exceptions_enable=False)

__version__ = "0.1.0"


@app.command()
def version():
    """
    Show the version and exit.
    """
    rprint(f"[bold green]RepeatsDB Lite 2, version: {__version__}[/bold green]\n")
    raise typer.Exit()


@app.command()
def query_directory(
        in_dir: Annotated[
            str, typer.Argument(help="Path to the input directory containing structure files (PDB/mmCIF)")
        ],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to the output directory")
        ],
        keep_temp: Annotated[bool, typer.Option(
            help="Keep temporary directory and files")] = False,
        max_eval: Annotated[
            float, typer.Option(help="Maximum E-value of the hits to prefilter")
        ] = cfg.max_eval,
        min_height: Annotated[
            float, typer.Option(
                help="Minimum height of TM-score signals to be processed")
        ] = cfg.min_height
):
    """
    Run the pipeline on a directory containing structure files (PDB/mmCIF).
    """
    temp_dir = os.path.join(out_dir, "temp")

    rprint(
        f"""[bold]Running RepeatsDB Lite 2 with the following parameters:[/bold]
[bold blue]Input directory: {in_dir} [/bold blue]
[bold blue]Output directory: {out_dir} [/bold blue]
[bold blue]Result directory: {temp_dir} [/bold blue]
[bold blue]Maximum E-value: {max_eval} [/bold blue]
[bold blue]Minimum height: {min_height} [/bold blue]\n"""
    )

    ex.execute_repeatsalgorithm(
        in_dir, out_dir, temp_dir, max_eval, min_height, keep_temp
    )


@app.command()
def download_pdb(
        pdb_id: Annotated[str, typer.Argument(help="PDB ID to download")],
        pdb_chain: Annotated[str, typer.Argument(help="PDB chain to query")],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to directory where output will be saved")
        ],
        keep_temp: Annotated[bool, typer.Option(
            help="Keep temporary directory and files")] = False,
        max_eval: Annotated[
            float, typer.Option(help="Maximum E-value of the hits to prefilter")
        ] = cfg.max_eval,
        min_height: Annotated[
            float, typer.Option(
                help="Minimum height of TM-score signals to be processed")
        ] = cfg.min_height
):
    """
    Download and query a structure from PDB by providing the PDB ID and the specific Chain of interest.
    """

    temp_dir = os.path.join(out_dir, "temp")

    in_dir = os.path.join(out_dir, "download_structures")

    success = ds.download_pdb_structure(pdb_id, pdb_chain, in_dir)
    if success == 1:
        raise typer.Abort()

    out_dir = os.path.join(out_dir, "out")

    rprint(
        f"""[bold]Running RepeatsDB Lite 2 with the following parameters:[/bold]
[bold blue]Input directory: {in_dir} [/bold blue]
[bold blue]Output directory: {out_dir} [/bold blue]
[bold blue]Result directory: {temp_dir} [/bold blue]
[bold blue]Maximum E-value: {max_eval} [/bold blue]
[bold blue]Minimum height: {min_height} [/bold blue]\n"""
    )

    ex.execute_repeatsalgorithm(
        in_dir, out_dir, temp_dir, max_eval, min_height, keep_temp
    )


@app.command()
def download_model(
        uniprot_id: Annotated[
            str, typer.Argument(
                help="UniProt ID of the AlphaFold structure to query")
        ],
        af_version: Annotated[
            str, typer.Argument(
                help="Version of AlphaFold to download structure from")
        ],
        out_dir: Annotated[
            str, typer.Argument(
                help="Path to directory where output will be saved")
        ],
        keep_temp: Annotated[bool, typer.Option(
            help="Keep temporary directory and files")] = False,
        max_eval: Annotated[
            float, typer.Option(help="Maximum E-value of the hits to prefilter")
        ] = cfg.max_eval,
        min_height: Annotated[
            float, typer.Option(
                help="Minimum height of TM-score signals to be processed")
        ] = cfg.min_height,
):
    """
    Download and query an AlphaFold model by providing the UniProt ID and the AlphaFold version of interest.
    """

    temp_dir = os.path.join(out_dir, "temp")

    in_dir = os.path.join(out_dir, "download_structures")

    rprint(f"""[bold]Downloading AlphaFold model for {uniprot_id}[/bold]\n""")

    success = ds.download_alphafold_structure(uniprot_id, af_version, in_dir)
    if success != 200:
        raise typer.Abort()

    out_dir = os.path.join(out_dir, "out")

    rprint(
        f"""[bold]Running RepeatsDB Lite 2 with the following parameters:[/bold]
[bold blue]Input directory: {in_dir} [/bold blue]
[bold blue]Output directory: {out_dir} [/bold blue]
[bold blue]Result directory: {temp_dir} [/bold blue]
[bold blue]Maximum E-value: {max_eval} [/bold blue]
[bold blue]Minimum height: {min_height} [/bold blue]\n"""
    )

    ex.execute_repeatsalgorithm(
        in_dir, out_dir, temp_dir, max_eval, min_height, keep_temp
    )


if __name__ == "__main__":
    app()
