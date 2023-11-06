import src.config as cfg
import os
import typer
from typing_extensions import Annotated
from src import execute_repeatsalgorithm as ex

# Initialize the app
app = typer.Typer()

__version__ = "0.1.0"


def version_callback(value: bool):
    if value:
        print(f"RepeatsDB Lite 2, version: {__version__}")
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
    ex.execute_repeatsalgorithm(
        in_dir, out_dir, result_dir, max_eval, min_height)


if __name__ == "__main__":
    app()
