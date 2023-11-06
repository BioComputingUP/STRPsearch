from . import search_vs_tul
import typer
import shutil
import mimetypes
import tempfile
import os
from . import config as cfg
from Bio.PDB import MMCIFParser, PDBIO, is_aa, PDBParser
import pandas as pd
import numpy as np
from .tmalign import Tmalign
from . import predictor_utils as utils
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

# Define the default hyperparameters
distance_p = cfg.distance_p
frame_step = cfg.frame_step
width_p = cfg.width_p
prominence_p = cfg.prominence_p
threshold_p = cfg.threshold_p
window_p = cfg.window_p
flexibility_p = 1 - distance_p

# Instantiate essential modules
cif_parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)
io_handler = PDBIO()
tmalign = Tmalign(exe_path=cfg.tmalign_exe_path)

# Specify the paths to RepeatsDB ground-truth library
tul_fs_db = "data/databases/tul_foldseek_db/db"
target_db = "data/databases/tmax"


def execute_repeatsalgorithm(
    structure_dir: str,
    out_dir: str,
    temp_dir: str,
    max_eval_p: float,
    height_p: float,
):
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        print("Output directory already exists")
        raise typer.Exit()
    try:
        os.mkdir(temp_dir)
    except FileExistsError:
        print("Results directory already exists")
        raise typer.Exit()
    fs_output = os.path.join(temp_dir, "fs_output.tsv")

    search_vs_tul.search_tul(
        foldseek_exe_path=cfg.foldseek_exe_path,
        query_dir=structure_dir,
        tul_fs_db=tul_fs_db,
        output_file=fs_output,
        temp_dir=temp_dir,
    )

    found_hit, target_df = search_vs_tul.find_target(
        output_file=fs_output, max_eval=max_eval_p
    )

    if found_hit:
        temp_query_dir_list = []
        for i in range(len(target_df)):
            row = target_df.iloc[i]
            # Define essential variables
            query_name = row["query"]
            target_name = row["target"]
            target_chain = target_name[4]
            e_val = row["e_value"]
            ct = str(row["t_ct"])
            target_avg_len = row["t_avg_length"]

            print("Query name:", query_name)
            print("Target name: ", target_name)
            print("Target CT: ", ct)

            insert_len = 60
            if ct == "4.4":
                insert_len *= 5

            query_path = os.path.join(os.path.join(structure_dir, query_name))
            mime_type, encoding = mimetypes.guess_type(query_path)
            if mime_type:
                if "pdb" in mime_type:
                    # Parse the structure from the PDB file
                    qstructure = pdb_parser.get_structure(
                        query_name, query_path)
                elif "cif" in mime_type:
                    # Parse the structure from the PDB file
                    qstructure = cif_parser.get_structure(
                        query_name, query_path)
                else:
                    print("Only PDB / mmCIF format is accepted for query files")
                    raise TypeError
            else:
                print("The query file format is ambiguous")
                raise TypeError

            query_chain = list(qstructure.get_chains())[0].id
            # Designate the chain from the structure
            qchain = qstructure[0][query_chain]
            # Extract residues from the chain
            qchain_residues = [
                residue
                for residue in qchain.get_residues()
                if is_aa(residue.get_resname())
            ]
            qchain_residue_nums = [chain_res.id[1]
                                   for chain_res in qchain_residues]

            target_tmax_path = os.path.join(
                target_db, target_name[1:3], f"{target_name}_tmax.pdb"
            )
            tstructure = pdb_parser.get_structure(
                target_name, target_tmax_path)
            tchain = tstructure[0][target_chain]

            total_fragment_list = utils.get_res_frames(
                qchain_residues, len(tchain), step=frame_step
            )

            fragment_dir = tempfile.mkdtemp()
            fragment_path_list = []
            for fragment in total_fragment_list:
                fragment_start = fragment[0].id[1]
                fragment_end = fragment[-1].id[1]
                fragment_out_path = os.path.join(
                    fragment_dir,
                    f"{query_name[:5]}_{fragment_start}_{fragment_end}.pdb",
                )
                utils.get_structure(
                    fragment, query_chain, qstructure, fragment_out_path, io_handler
                )
                fragment_path_list.append(fragment_out_path)

            tmscore_results = []
            for fragment_path in fragment_path_list:
                fragment_start = os.path.basename(fragment_path).split("_")[1]
                command_output = tmalign(target_tmax_path, fragment_path)
                lines = command_output.split("\n")
                tm_score = None
                for line in lines:
                    if (
                        line.startswith("TM-score=")
                        and "if normalized by average length" in line
                    ):
                        # Extract the TM-score value
                        tm_score = float(line.split("=")[1].split()[0])
                        break  # Stop looping after finding the desired TM-score line

                tmscore_results.append([fragment_start, tm_score])

            tm_score_results = [[int(i[0]), i[1]] for i in tmscore_results]
            tm_score_results.sort(key=lambda x: x[0])

            x = [i[0] for i in tm_score_results]
            y = [i[1] for i in tm_score_results]

            y = np.array(y)
            x = np.array(x)

            temp_query_dir = os.path.join(temp_dir, query_name)
            os.makedirs(temp_query_dir, exist_ok=True)
            temp_query_dir_list.append(temp_query_dir)

            window_size = round(target_avg_len * window_p)
            if window_size % 2 == 0:
                window_size -= 1

            poly_order = 0
            y = savgol_filter(y, window_size, poly_order)

            x, y = utils.adjust_graph_ends(x=x, y=y, frame_step=frame_step)

            peaks, _ = find_peaks(
                y,
                height=height_p,
                distance=round(len(tchain) * distance_p),
                width=width_p,
                prominence=prominence_p,
                threshold=threshold_p,
            )
            start_list = [
                x[peak_idx] for peak_idx in peaks if x[peak_idx] in qchain_residue_nums
            ]

            predicted_unit_list = utils.calculate_ranges(
                start_list,
                round(target_avg_len),
                qchain_residues[-1].id[1],
                flexibility=flexibility_p,
            )

            regions_dict = utils.parse_regions(
                query_name, predicted_unit_list, insert_len
            )

            if regions_dict:
                out_name = f"{query_name}_{target_name}_{ct}_{e_val}"

                figure_outpath = os.path.join(
                    temp_query_dir, f"{out_name}.png")
                utils.plot_tmscore_graph(
                    x=x, y=y, regions_dict=regions_dict, out_path=figure_outpath
                )

                json_outpath = os.path.join(temp_query_dir, f"{out_name}.json")
                utils.make_json(
                    pdb_id_chain=query_name,
                    ct=ct,
                    regions_dict=regions_dict,
                    out_path=json_outpath,
                )

                for region_id, v in regions_dict.items():
                    start_res = v["units"][0][0]
                    end_res = v["units"][-1][1]

                    region_outpath = os.path.join(
                        temp_query_dir, f"{region_id}_{target_name}_{ct}_{e_val}.pdb"
                    )
                    region_range = utils.get_chain_range(
                        start_res, end_res, qchain_residues
                    )
                    utils.get_structure(
                        region_range,
                        query_chain,
                        qstructure,
                        region_outpath,
                        io_handler,
                    )

                shutil.rmtree(fragment_dir)
                # shutil.rmtree(target_dir)
            else:
                print(
                    f"No repeat region was found with the spcified min_height of {height_p}"
                )
        for temp_query_dir in temp_query_dir_list:
            dir_name = os.path.basename(temp_query_dir)
            out_query_dir = os.path.join(out_dir, dir_name)
            os.makedirs(out_query_dir, exist_ok=True)

            filename_dict = {"q_start": [], "q_end": [], "e_value": []}
            for filename in os.listdir(temp_query_dir):
                if filename.endswith(".pdb"):
                    print(filename)
                    filename_motifs = filename[len(query_name):].split("_")
                    q_start = int(filename_motifs[1])
                    q_end = int(filename_motifs[2])
                    e_value = float(filename_motifs[-1].strip(".pdb"))
                    filename_dict["q_start"].append(q_start)
                    filename_dict["q_end"].append(q_end)
                    filename_dict["e_value"].append(e_value)

            filename_df = pd.DataFrame(filename_dict)
            filtered_df = search_vs_tul.filter_overlap(filename_df)

            src_filepaths = []
            for filename in os.listdir(temp_query_dir):
                for e_value in filtered_df["e_value"]:
                    if str(e_value) in filename:
                        filepath = os.path.join(temp_query_dir, filename)
                        src_filepaths.append(filepath)

            for filepath in src_filepaths:
                filename = os.path.basename(filepath)
                dst_path = os.path.join(out_query_dir, filename)
                shutil.copy(filepath, dst_path)

        shutil.rmtree(temp_dir)
    else:
        print(f"No hit was found below the specified max_eval of {max_eval_p}")
