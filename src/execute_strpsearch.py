import json
from Bio.PDB import MMCIFParser, PDBIO, is_aa, PDBParser
from scipy.signal import find_peaks
from rich import print as rprint
from . import alignment_utils as au
from . import general_utils as gu
from . import config as cfg
import pandas as pd
import traceback
import tempfile
import logging
import shutil
import typer
import os
import time
import re
import gemmi

# Define computational parameters
distance_p = cfg.distance_p
frame_step = cfg.frame_step_p
width_p = cfg.width_p
prominence_p = cfg.prominence_p
threshold_p = cfg.threshold_p
window_p = cfg.window_p
flexibility_p = 1 - distance_p

# Instantiate essential modules
cif_parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)
io_handler = PDBIO()

# Specify paths to ground-truth libraries
tul_db = "data/databases/tul_foldseek_db/db"
rul_db = "data/databases/rul_structure_db/"
ontology_df = pd.read_csv("data/ontology.tsv", delimiter="\t")

# Load classification TM-score data
with open("data/ct_tmscore_means.json", 'r') as fp:
    ct_tmscore_dict = json.load(fp)


def execute_predstrp(
        structure_dir: str,
        out_dir: str,
        temp_dir: str,
        keep_temp: bool,
        pymol_pse: bool,
        max_eval_p: float,
        min_height_p: str,
):
    """
    Executes the PredSTRP pipeline to identify and analyze repeat regions in protein structures.

    Args:
        structure_dir (str): Directory containing input structures.
        out_dir (str): Directory to save output results.
        temp_dir (str): Temporary directory for intermediate files.
        keep_temp (bool): Whether to keep temporary files.
        pymol_pse (bool): Whether to generate PyMOL session files.
        max_eval_p (float): Maximum E-value threshold for filtering hits.
        min_height_p (str): Minimum height for peak detection ("auto" or a float value).

    Returns:
        None
    """
    start_time = time.time()

    # Configure logging
    logging.basicConfig(
        filename=os.path.join(out_dir, "debug.log"),
        filemode="w",
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        level=logging.WARNING
    )

    # Validate input files
    pdb_cif_fps, have_space_fps = gu.check_files(in_dir=structure_dir)
    if len(pdb_cif_fps) > 0:
        if len(have_space_fps) > 0:
            rprint(f"\n[bold][{gu.time()}][/bold] [bold red]"
                   f"Please remove spaces from the following files:\n")
            logging.error("Files with spaces detected:")
            for fp in have_space_fps:
                rprint(f"[bold red]{fp}\n")
                logging.error(fp)
            raise typer.Abort()
    else:
        rprint(f"\n[bold][{gu.time()}][/bold] [bold red]"
               f"No file with pdb/cif extension was found")
        logging.error("No file with pdb/cif extension found")
        raise typer.Abort()

    # Display pipeline parameters
    rprint(f"\n[bold][{gu.time()}] Running PredSTRP with the following parameters:")
    rprint(f"[bold blue]Input directory: {structure_dir}")
    rprint(f"[bold blue]Output directory: {out_dir}")
    rprint(f"[bold blue]Temporary directory: {temp_dir}")
    rprint(f"[bold blue]Maximum E-value: {max_eval_p}")
    rprint(f"[bold blue]Minimum height: {min_height_p}\n")

    # Specify the path to save Foldseek search output
    fs_output = os.path.join(temp_dir, "fs_output.tsv")

    # Perform Foldseek search
    rprint(f"[bold][{gu.time()}] Finding potential hits ...\n")
    au.search_tul(
        foldseek_exe_path=cfg.foldseek_exe_path,
        query_dir=structure_dir,
        tul_fs_db=tul_db,
        output_file=fs_output,
        temp_dir=temp_dir,
    )

    # Parse Foldseek results
    found_hit, target_df = au.find_target(output_file=fs_output, max_eval=max_eval_p)
    print("The target hits obtained from Foldseek:\n", target_df)

    # Process hits if found
    if found_hit:
        temp_query_dir_list = []
        error_count = 0

        for idx in range(len(target_df)):
            try:
                row = target_df.iloc[idx]
                query_id = "_".join(row["query"].split("_")[:-1])
                query_chain = re.search(r"_([^_]+)$", row["query"]).group(1)
                query_name = query_id + "_" + query_chain
                target_name = row["target"]
                target_chain = target_name[4]
                e_val = row["e_value"]
                ct = str(row["t_ct"])
                target_avg_len = float(row["t_avg_length"])
                target_classi = gu.get_repeat_classi(ontology_df=ontology_df, code=ct)

                max_insertion_p = 60
                if ct == "4.4":
                    max_insertion_p *= 2

                # Determine minimum height for peak detection
                if min_height_p == "auto":
                    if ct in ct_tmscore_dict.keys():
                        final_min_height_p = round(ct_tmscore_dict[ct][0] - ct_tmscore_dict[ct][1], 2)
                    else:
                        final_min_height_p = 0.35
                else:
                    final_min_height_p = float(min_height_p)

                rprint(f"[bold][{gu.time()}] Processing hit {idx + 1}/{len(target_df)} ...")
                rprint(f"[bold blue]Query: {query_name}")
                rprint(f"[bold blue]Target: {target_name}")
                rprint(f"[bold blue]Classification: {target_classi}\n")

                # Load query structure and chain
                query_path = os.path.join(structure_dir, f"{query_name}.cif")
                qstructure = gemmi.read_structure(query_path)
                qmodel = qstructure[0]
                qchain_letter = qmodel[0].name
                qchain = qmodel[qchain_letter]
                qchain_residues = [res for res in qchain if gemmi.find_tabulated_residue(res.name)]
                qchain_residue_nums = [res.seqid.num for res in qchain_residues]

                # Load target representative unit
                target_repunit_path = os.path.join(rul_db, target_name[1:3], f"{target_name}_tmax.pdb")
                tstructure = pdb_parser.get_structure(target_name, target_repunit_path)
                tchain = tstructure[0][target_chain]

                # Fragmentize query based on target unit length
                total_fragment_list = gu.get_res_frames(
                    res_range=qchain_residues,
                    length=len(tchain),
                    step=frame_step
                )

                # Save fragments to temporary directory
                fragment_dir = tempfile.mkdtemp(dir=temp_dir)
                fragment_path_list = []
                for fragment in total_fragment_list:
                    fragment_start = fragment[0].seqid.num
                    fragment_end = fragment[-1].seqid.num
                    fragment_out_path = os.path.join(fragment_dir, f"{query_name}_{fragment_start}_{fragment_end}.cif")
                    gu.get_structure_cif(
                        res_range=[res.seqid.num for res in fragment],
                        chain_id=qchain_letter,
                        structure=qstructure,
                        out_path=fragment_out_path
                    )
                    fragment_path_list.append(fragment_out_path)

                # Generate TM-score graph data
                x, y = au.get_tmscore_graph_data_us(
                    query_name=query_name,
                    fragment_path_list=fragment_path_list,
                    target_repunit_path=target_repunit_path,
                    usalign_exe_path=cfg.usalign_exe_path
                )

                # Smooth and adjust graph data
                x, y = list(x), list(y)
                y = gu.smooth_graph(y=y, target_avg_len=target_avg_len, window_p=window_p)
                x, y = gu.adjust_graph_ends(x=x, y=y, frame_step=frame_step)

                # Detect peaks in the graph
                peaks, _ = find_peaks(
                    x=y,
                    height=final_min_height_p,
                    distance=round(len(tchain) * distance_p),
                    width=width_p,
                    prominence=prominence_p,
                    threshold=threshold_p,
                )
                peak_residue_nums = [x[peak_idx] for peak_idx in peaks if x[peak_idx] in qchain_residue_nums]

                # Calculate repeat regions based on peaks
                predicted_unit_list = gu.calculate_ranges(
                    peak_residue_nums=peak_residue_nums,
                    distance=round(target_avg_len),
                    max_length=qchain_residues[-1].seqid.num,
                    flexibility=flexibility_p
                )

                # Map and save repeat regions
                regions_dict = gu.parse_regions(query_name, predicted_unit_list, max_insertion_p)
                if regions_dict:
                    rprint(f"[bold][{gu.time()}][/bold] [bold green]"
                           f"{len(regions_dict)} potential repeat regions found with hit {idx + 1}/{len(target_df)}")
                    # Save regions and generate outputs
                    # (Details omitted for brevity)
                else:
                    rprint(f"[bold][{gu.time()}][/bold] [bold yellow]"
                           f"No potential repeat region found with hit {idx + 1}/{len(target_df)}\n")

                # Clean up temporary fragment directory
                shutil.rmtree(fragment_dir)

            except Exception as e:
                error_count += 1
                traceback.print_exc()
                logging.error(traceback.format_exc())

        # Finalize results and clean up
        # (Details omitted for brevity)

    else:
        # Handle case where no hits were found
        rprint(f"[bold][{gu.time()}][/bold] [bold yellow]"
               f"No hit found below the specified maximum E-value of {max_eval_p}")
        shutil.rmtree(temp_dir)
