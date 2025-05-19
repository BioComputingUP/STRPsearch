from src.foldseek import Foldseek
import pandas as pd
import numpy as np
import os
from src.usalign import Usalign
from Bio.PDB import Select




def search_tul(foldseek_exe_path, query_dir, tul_fs_db, output_file, temp_dir):
    """
    Searches the PDB/mmCIF structures in a directory against a Foldseek database and saves the output.
    """
    foldseek = Foldseek(exe_path=foldseek_exe_path)
    columns = "query,target,evalue,qstart,qend"
    foldseek.easy_search(
        input_file=query_dir,
        db=tul_fs_db,
        output_file=output_file,
        temp_dir=temp_dir,
        columns=columns,
    )


def check_overlap(range1, range2):
    """
    Checks if two ranges have at least 50% overlap and returns a boolean accordingly.
    """
    range1 = set(range1)
    overlapped = len(range1.intersection(range2))
    overlap_percentage = overlapped / len(range1)
    return overlap_percentage >= 0.5


def filter_overlap(df):
    """
    Filters a dataframe to retain entries with the smallest E-value when ranges overlap by more than 50%.
    """
    df_pieces = []
    for query in df["query"].unique():
        query_df = df[df["query"] == query].reset_index()
        to_keep_indices = [query_df["e_value"].idxmin()]

        for q_i in range(len(query_df)):
            if q_i not in to_keep_indices:
                q_row = query_df.iloc[q_i]
                q_eval = q_row["e_value"]
                q_range = range(int(q_row["q_start"]), int(q_row["q_end"]))
                overlap = []

                for idx, master_i in enumerate(to_keep_indices):
                    master_row = query_df.iloc[master_i]
                    master_eval = master_row["e_value"]
                    master_range = range(int(master_row["q_start"]), int(master_row["q_end"]))

                    if check_overlap(master_range, q_range):
                        overlap.append(True)
                        if q_eval < master_eval:
                            to_keep_indices[idx] = q_i
                    else:
                        overlap.append(False)

                if True not in overlap:
                    to_keep_indices.append(q_i)

        df_piece = query_df.iloc[to_keep_indices]
        df_pieces.append(df_piece)

    return pd.concat(df_pieces, ignore_index=True).sort_values(by="query")


def find_target(output_file, max_eval):
    """
    Parses and filters the Foldseek search output to retain hits with the highest potential.
    """
    df = pd.read_csv(output_file, delimiter="\t", header=None, dtype={"ctr": str})
    df.columns = ["query", "target", "e_value", "q_start", "q_end"]
    df = df[df["e_value"] <= max_eval]

    if len(df) == 0:
        return False, df

    # Extract class topology (ct) and average length from the target column
    df[['t_ct', 't_avg_length']] = df['target'].str.extract(r'_(\d+\.\d+)_(\d+\.\d+)\.pdb')
    df['t_ct'] = df['t_ct'].astype(str)
    df['t_avg_length'] = df['t_avg_length'].astype(float)
    df['target'] = df['target'].apply(lambda x: '_'.join(x.split('_')[:-2]))

    target_list = []

    for query in df["query"].unique():
        query_df = df[df["query"] == query]
        lowest_eval_indices = query_df.groupby('t_ct')["e_value"].idxmin()
        lowest_eval_rows = query_df.loc[lowest_eval_indices]

        for _, row in lowest_eval_rows.iterrows():
            target_list.append({
                "query": row["query"],
                "target": row["target"],
                "q_start": row["q_start"],
                "q_end": row["q_end"],
                "t_ct": row["t_ct"],
                "t_avg_length": row["t_avg_length"],
                "e_value": row["e_value"],
            })

    target_df = pd.DataFrame(target_list)
    return True, filter_overlap(target_df)


def get_tmscore_graph_data(query_name, fragment_path_list, target_repunit_path, tmalign_exe_path):
    """
    Generates TM-score graph data by aligning structure fragments with the target structure using TM-align.
    """
    tmalign = Tmalign(exe_path=tmalign_exe_path)
    tmscore_results = []

    for fragment_path in fragment_path_list:
        filename = os.path.basename(fragment_path)
        fragment_start = filename[len(query_name):].split("_")[1]
        command_output = tmalign(target_repunit_path, fragment_path)

        # Parse TM-score from TM-align output
        tm_score = None
        for line in command_output.split("\n"):
            if line.startswith('TM-score=') and 'if normalized by average length' in line:
                tm_score = float(line.split('=')[1].split()[0])
                break

        tmscore_results.append([fragment_start, tm_score])

    # Sort results by fragment start position
    tm_score_results = [[int(i[0]), i[1]] for i in tmscore_results]
    tm_score_results.sort(key=lambda x: x[0])

    x = np.array([i[0] for i in tm_score_results])
    y = np.array([i[1] for i in tm_score_results])

    return x, y


def get_tmscore_graph_data_us(query_name, fragment_path_list, target_repunit_path, usalign_exe_path):
    """
    Generates TM-score graph data by aligning structure fragments with the target structure using US-align.
    """
    usalign = Usalign(exe_path=usalign_exe_path)
    tmscore_results = []

    for fragment_path in fragment_path_list:
        filename = os.path.basename(fragment_path)
        fragment_start = filename[len(query_name):].split("_")[1]
        command_output = usalign(target_repunit_path, fragment_path)

        # Parse TM-score from US-align output
        tm_score = None
        for line in command_output.splitlines():
            if line.startswith('TM-score='):
                tm_score = float(line.split('=')[1].split()[0])
                break

        tmscore_results.append([fragment_start, tm_score])

    # Sort results by fragment start position
    tm_score_results = [[int(i[0]), i[1]] for i in tmscore_results]
    tm_score_results.sort(key=lambda x: x[0])

    x = np.array([i[0] for i in tm_score_results])
    y = np.array([i[1] for i in tm_score_results])

    return x, y