import json
from Bio.PDB import MMCIFParser, PDBIO, MMCIFIO, is_aa, PDBParser
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
import warnings
from Bio import BiopythonWarning
import time
import re
import gemmi

warnings.filterwarnings("ignore", category=BiopythonWarning)
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
io_handler_cif = MMCIFIO()

project_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
ontology_path = os.path.join(project_root, "data", "ontology.tsv")
ontology_df = pd.read_csv(ontology_path, delimiter="\t")

ct_tmscore_path = os.path.join(project_root, "data", "ct_tmscore_means.json")
# Load classification TM-score data
with open(ct_tmscore_path, 'r') as fp:
    ct_tmscore_dict = json.load(fp)


def execute_predstrp(
        structure_dir: str,
        out_dir: str,
        temp_dir: str,
        keep_temp: bool,
        pymol_pse: bool,
        max_eval_p: float,
        min_height_p: str,
        tul_db: str ,
        rul_db: str,
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

    
    gu.segment_cif_directory(structure_dir,structure_dir)
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
    # Process hits if found
    if found_hit:
        temp_query_dir_list = []
        error_count = 0

        for idx in range(len(target_df)):
            try:
                row = target_df.iloc[idx]
                query_name_path= row["query"]
                parts= query_name_path.split("_")
                query_id =parts[0]
                query_chain = parts[1]
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
                query_path = os.path.join(structure_dir, f"{query_name_path}.cif")
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
                if x is None or y is None or len(x) == 0 or len(y) == 0:
                    rprint(f"[bold][{gu.time()}][/bold] [bold yellow]"
                           f"TM-score graph data is empty for hit {idx + 1}/{len(target_df)}. Skipping.\n")
                    continue
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
                    if len(regions_dict) > 1:
                        rprint(f"[bold][{gu.time()}][/bold] [bold green]"
                               f"{len(regions_dict)} potential repeat regions were found with hit "
                               f"{idx + 1}/{len(target_df)}:")
                    else:
                        rprint(f"[bold][{gu.time()}][/bold] [bold green]"
                               f"1 potential repeat region was found with hit "
                               f"{idx + 1}/{len(target_df)}:")

                    for idx, region_id in enumerate(list(regions_dict.keys())):
                        rprint(f"[bold blue]"
                               f"Region {idx + 1}: {region_id}")

                    # Create a directory with the name of the query at the designated temporary directory
                    temp_query_dir = os.path.join(temp_dir, f"{query_name}")
                    os.makedirs(temp_query_dir, exist_ok=True)
                    temp_query_dir_list.append(temp_query_dir)

                    # Loop the through the regions
                    for region_id, components in regions_dict.items():
                        start_res = components["units"][0][0]
                        end_res = components["units"][-1][1]

                        # Because the end residue of each repeat region is always estimated
                        # Check if it actually exists on the query chain and not missing
                        # If missing (index error), find the largest smaller residue number before (error handling)
                        try:
                            gu.get_res_index(end_res, qchain_residues)
                        except IndexError:
                            end_res = gu.find_largest_smaller_number(end_res, qchain_residues)

                        # Define an output name
                        out_name = f"{region_id}_{ct}_{e_val}"

                        region_out_path = os.path.join(temp_query_dir, f"{out_name}.cif")

                        # Extract and save the structure of the region
                        region_range = gu.get_chain_range(
                            start=start_res,
                            end=end_res,
                            chain_residues=qchain_residues
                        )
                        gu.get_structure_cif(
                            res_range=region_range,
                            chain_id=qchain_letter,
                            structure=qstructure,
                            out_path=region_out_path,
                        )

                        # Create and save the PyMOL session of the repeat region highlighted by the integral components
                        if pymol_pse:
                            pymol_out_path = os.path.join(temp_query_dir, f"{out_name}.pse")

                            gu.create_pymol_session(
                                region_id=region_id,
                                structure_path=query_path,
                                components=components,
                                output_path=pymol_out_path
                            )

                        # Create the mapped graph plot
                        figure_out_path = os.path.join(temp_query_dir, f"{out_name}.png")

                        gu.plot_tmscore_graph(
                            x=x,
                            y=y,
                            region_components=components,
                            out_path=figure_out_path
                        )

                        # Create the JSON file
                        json_out_path = os.path.join(temp_query_dir, f"{out_name}.json")

                        gu.make_json(
                            structure_id=query_name,
                            chain_id=qchain_letter,
                            ct=ct,
                            region_id=region_id,
                            components=components,
                            out_path=json_out_path,
                        )

                    # Delete the fragment directory
                    shutil.rmtree(fragment_dir)

                else:
                    rprint(f"[bold][{gu.time()}][/bold] [bold yellow]"
                           f"No potential repeat region was found with hit "
                           f"{idx + 1}/{len(target_df)}\n")
            except Exception as e:
                error_count += 1
                # traceback.print_exc()
                logging.error(traceback.format_exc())
                rprint(f"[bold yellow]WARNING: Error occurred in hit {idx+1}/{len(target_df)}: {e}[/bold yellow]")

        if temp_query_dir_list:
            rprint(f"\n[bold][{gu.time()}] "
                   f"All the hits were processed\n")

            # Transfer the files from temporary directory to the final query output directory
            rprint(f"[bold]"
                   f"[{gu.time()}] "
                   f"Transferring the final results to the output directory ...\n")
        

        # Loop the output directories in the temporary directory
        for temp_query_dir in temp_query_dir_list:
            try:
                
                query_name = "_".join(os.path.basename(temp_query_dir).split("_")[:2])
                # Create a dict to save certain attributes of the mapped regions
                # Loop through the outputs in the temporary directory and extract attributes
                filename_dict = {"query": [], "q_start": [], "q_end": [], "e_value": []}
                for filename in os.listdir(temp_query_dir):
                    if filename.endswith(".cif"):
                        filename_motifs = filename[len(query_name):].split("_")
                        q_start = int(filename_motifs[1])
                        q_end = int(filename_motifs[2])
                        e_value = float(filename_motifs[-1].strip(".cif"))
                        filename_dict["query"].append(query_name)
                        filename_dict["q_start"].append(q_start)
                        filename_dict["q_end"].append(q_end)
                        filename_dict["e_value"].append(e_value)

                filename_df = pd.DataFrame(filename_dict)
                # Filter the overlapping regions
                filtered_df = au.filter_overlap(filename_df)

                src_region_fps_dict = {}
                for row_idx in range(len(filtered_df)):
                    region_num = str(row_idx + 1)
                    row = filtered_df.iloc[row_idx]
                    query_name = row["query"]
                    q_start = row["q_start"]
                    q_end = row["q_end"]
                    basename = f"{query_name}_{q_start}_{q_end}"
                    for filename in os.listdir(temp_query_dir):
                        if basename in filename:
                            filepath = os.path.join(temp_query_dir, filename)
                            if region_num in src_region_fps_dict.keys():
                                src_region_fps_dict[region_num].append(filepath)
                            else:
                                src_region_fps_dict[region_num] = [filepath]

                # Get the name of the current query directory
                dir_name = os.path.basename(temp_query_dir)
                # For each query, define the path of the final output directory of that query
                out_query_dir = os.path.join(
                    out_dir,
                    f"{'_'.join(dir_name.split('_')[:-1])}_results",
                    f"chain_{dir_name.split('_')[-1]}"
                )
                # Create the final query output directory
                os.makedirs(out_query_dir, exist_ok=True)
                # Convert the dict to pandas df

                for region_num, filepaths in src_region_fps_dict.items():
                    out_region_dir = os.path.join(out_query_dir, f"region_{region_num}")
                    os.makedirs(out_region_dir, exist_ok=True)
                    for filepath in filepaths:
                        extension = os.path.basename(filepath).split(".")[-1]
                        dst_name = "_".join(os.path.basename(filepath).split("_")[:-2]) + "." + extension
                        dst_path = os.path.join(out_region_dir, dst_name)
                        shutil.copy(filepath, dst_path)

            except Exception as e:
                error_count += 1
                # traceback.print_exc()
                logging.error(traceback.format_exc())
                rprint(f"[bold yellow]WARNING: Error occurred while transferring files for query directory '{temp_query_dir}': {e}[/bold yellow]")

        end_time = time.time()
        elapsed_time = end_time - start_time
        time_stamp_path = os.path.join(out_dir, "time.txt")
        with open(time_stamp_path, "w") as fp:
            fp.write(str(elapsed_time))

        rprint(f"[bold][{gu.time()}] "
               f"Task finished with {error_count} errors")

        # If keep_temp file is set to False, delete the temporary directory
        if not keep_temp:
            shutil.rmtree(temp_dir)

    # If no hit was found, remove the created temp and output directories
    else:
        end_time = time.time()
        elapsed_time = end_time - start_time
        time_stamp_path = os.path.join(out_dir, "time.txt")
        with open(time_stamp_path, "w") as fp:
            fp.write(str(elapsed_time))
        rprint(f"[bold][{gu.time()}][/bold] [bold yellow]"
               f"No hit was found below the specified maximum E-value of {max_eval_p}[/bold yellow]")
        shutil.rmtree(temp_dir)
