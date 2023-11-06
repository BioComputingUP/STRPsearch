import src.config as cfg
import os
import typer
from typing_extensions import Annotated
from src import execute_repeatsalgorithm as ex
from rich import print
from src import download_structure as ds

# Initialize the app
app = typer.Typer()

__version__ = "0.1.0"


def version_callback(value: bool):
    if value:
        print(
            f"[bold green]RepeatsDB Lite 2, version: {__version__}[/bold green]")
        raise typer.Exit()


@app.command()
def directory(
    in_dir: Annotated[
        str, typer.Argument(help="Path to directory containing PDB files")
    ],
    out_dir: Annotated[
        str, typer.Argument(
            help="Path to directory where output will be saved")
    ],
    result_dir: Annotated[
        str,
        typer.Argument(
            help="Path to directory where result files will be saved. Defaults to a new directory named 'results' in the output directory"
        ),
    ] = "results",
    max_eval: Annotated[
        float, typer.Option(help="Maximum E-value of the targets to prefilter")
    ] = cfg.max_eval,
    min_height: Annotated[
        float, typer.Option(
            help="Minimum height of TM-score signals to be processed")
    ] = cfg.min_height,
    version: Annotated[
        bool,
        typer.Option(help="Show tool version", callback=version_callback),
    ] = False,
):
    """
    Run the pipeline on a directory containing PDB files.
    """
    result_dir = os.path.join(out_dir, result_dir)

    print(
        f"""\n[bold]Running RepeatsDB Lite 2 with the following parameters:[/bold]
[bold blue]Input directory: {in_dir} [/bold blue]
[bold blue]Output directory: {out_dir} [/bold blue]
[bold blue]Result directory: {result_dir} [/bold blue]
[bold blue]Maximum E-value: {max_eval} [/bold blue]
[bold blue]Minimum height: {min_height} [/bold blue]\n"""
    )

    ex.execute_repeatsalgorithm(
        in_dir, out_dir, result_dir, max_eval, min_height)


@app.command()
def download_pdb(
    pdb_id: Annotated[str, typer.Argument(help="PDB ID to download")],
    pdb_chain: Annotated[str, typer.Argument(help="PDB chain to query")],
    out_dir: Annotated[
        str, typer.Argument(
            help="Path to directory where output will be saved")
    ],
    result_dir: Annotated[
        str,
        typer.Argument(
            help="Path to directory where result files will be saved. Defaults to a new directory named 'results' in the output directory"
        ),
    ] = "results",
    max_eval: Annotated[
        float, typer.Option(help="Maximum E-value of the targets to prefilter")
    ] = cfg.max_eval,
    min_height: Annotated[
        float, typer.Option(
            help="Minimum height of TM-score signals to be processed")
    ] = cfg.min_height,
    version: Annotated[
        bool,
        typer.Option(help="Show tool version", callback=version_callback),
    ] = False,
):
    """
    Run the pipeline downloading a structure and querying a specific chain.
    """

    result_dir = os.path.join(out_dir, result_dir)

    in_dir = os.path.join(out_dir, "download_structures")

    print(f"""\n[bold]Downloading PDB structure {pdb_id}[/bold]""")

    ds.download_pdb_structure(pdb_id, pdb_chain, in_dir)

    out_dir = os.path.join(out_dir, "out")

    print(
        f"""\n[bold]Running RepeatsDB Lite 2 with the following parameters:[/bold]
[bold blue]Input directory: {in_dir} [/bold blue]
[bold blue]Output directory: {out_dir} [/bold blue]
[bold blue]Result directory: {result_dir} [/bold blue]
[bold blue]Maximum E-value: {max_eval} [/bold blue]
[bold blue]Minimum height: {min_height} [/bold blue]\n"""
    )

    ex.execute_repeatsalgorithm(
        in_dir, out_dir, result_dir, max_eval, min_height)


if __name__ == "__main__":
    app()
