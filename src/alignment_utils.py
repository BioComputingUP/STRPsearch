from src.foldseek import Foldseek
import pandas as pd
import numpy as np
import os
from src.tmalign import Tmalign


def search_tul(foldseek_exe_path, query_dir, tul_fs_db, output_file, temp_dir):
    """
    Searches the PDB/mmCIF structures in a directory against a Foldseek database and saves the output
    """

    # Initiate Foldseek
    foldseek = Foldseek(exe_path=foldseek_exe_path)
    # Define the desired columns to obtain from search
    columns = "query,target,evalue,qstart,qend"
    # Do the search by easy_search module
    foldseek.easy_search(
        input_file=query_dir,
        db=tul_fs_db,
        output_file=output_file,
        temp_dir=temp_dir,
        columns=columns,
    )


def check_overlap(range1, range2):
    """
    Checks if two ranges have at least 50% overlap and returns a BOOL accordingly
    """

    range1 = set(range1)
    overlapped = len(range1.intersection(range2))
    overlap_percentage = overlapped / len(range1)
    if overlap_percentage >= 0.66:
        return True
    return False


def filter_overlap(df):
    """
    Takes a dataframe, compares the entries 1-vs-1, if they overlap more than 50% -->
    keep the one that has the smallest E-value
    """

    # Create a list to save pieces of the df that have been filtered individually
    df_pieces = []
    # Loop the unique queries in the df
    for query in df["query"].unique():
        # Isolate the df based on the query
        query_df = df[df["query"] == query].reset_index()
        # Create a list to save the query_df row indices that want to keep
        to_keep_indices = []
        # Add the index of the row with the least E-value to the to_keep_indices
        min_eval_index = query_df["e_value"].idxmin()
        to_keep_indices.append(min_eval_index)
        # Loop the range of the df length (row indices)
        for q_i in range(len(query_df)):
            # Ignore the index if already in to_keep_indic
            if q_i not in to_keep_indices:
                # Access the row and extract the essential data regarding the range and E-value
                q_row = query_df.iloc[q_i]
                q_eval = q_row["e_value"]
                q_range = range(int(q_row["q_start"]), int(q_row["q_end"]))
                # Create a BOOL list to check if the row range overlaps with any from to_keep_indices
                overlap = []
                # Compare the range and E-value to the ones in to_keep_indices
                for idx, master_i in enumerate(to_keep_indices):
                    master_row = query_df.iloc[master_i]
                    master_eval = master_row["e_value"]
                    master_range = range(int(master_row["q_start"]), int(master_row["q_end"]))

                    # Check if they overlap at least 50%
                    if check_overlap(master_range, q_range):
                        # If yes, append True to overlap list, and if the current one has a lower E-value, replace
                        overlap.append(True)
                        if q_eval < master_eval:
                            to_keep_indices[idx] = q_i
                    # If not, ignore it and append False to the overlap list
                    else:
                        overlap.append(False)
                # If no overlap after all, added it to the to_keep_indices
                if True not in overlap:
                    to_keep_indices.append(q_i)
        # Extract only the rows associated with the kept indices in the query df
        df_piece = query_df.iloc[to_keep_indices]
        # Save it as a piece to df_pieces list
        df_pieces.append(df_piece)
        # Concat the pieces, sort based on the query names, and return it
    filtered_sorted_df = pd.concat(df_pieces, ignore_index=True).sort_values(by="query")
    return filtered_sorted_df


def find_target(output_file, max_eval):
    """
    Takes the TSV output file of Foldseek search,
    parses and filters the df to keep the hits with the highest potential
    """

    # Load the dataframe
    df = pd.read_csv(output_file, delimiter="\t", header=None, dtype={"ctr": str})
    # Name the columns
    df.columns = ["query", "target", "e_value", "q_start", "q_end"]
    # Filter the column based on the Maximum E-value allowed
    df = df[df["e_value"] <= max_eval]

    # Return False in case of no hits
    if len(df) == 0:
        return False, df

    # Extracting 'avg_len' and 'ct' from 'target' and assign them to new columns
    df[['t_ct', 't_avg_length']] = df['target'].str.extract(r'_(\d+\.\d+)_(\d+\.\d+)\.pdb')

    # Convert 'ct' and 'avg_len' to floating-point numbers
    df['t_ct'] = df['t_ct'].astype(str)
    df['t_avg_length'] = df['t_avg_length'].astype(float)

    # Modify 'target' column
    df['target'] = df['target'].apply(lambda x: '_'.join(x.split('_')[:-2]))

    # Create a list to save the query-target combos that have the least E-value for each class.topology (ct)
    target_list = []
    for query in df["query"].unique():
        query_df = df[df["query"] == query]
        # Group the DataFrame by unique class topology values and find the index with the lowest E-value for each group
        lowest_eval_indices = query_df.groupby('t_ct')["e_value"].idxmin()

        # Extract the rows with the lowest eval for each class topology
        lowest_eval_rows = query_df.loc[lowest_eval_indices]
        for i in lowest_eval_rows.index:
            row = lowest_eval_rows.loc[i]
            target_list.append({"query": row["query"], "target": row["target"], "q_start": row["q_start"],
                                "q_end": row["q_end"], "t_ct": row["t_ct"], "t_avg_length": row["t_avg_length"],
                                "e_value": row["e_value"]})

    # Convert the list to a dataframe
    target_df = pd.DataFrame(target_list)
    # Subject it to the filter_overlap function
    target_df = filter_overlap(target_df)

    return True, target_df


def get_tmscore_graph_data(query_name, fragment_path_list, target_repunit_path, tmalign_exe_path):
    """
    Generates the TM-score graph data by aligning consequent structure fragments with the target structure
    """

    # Initiate TM-align
    tmalign = Tmalign(exe_path=tmalign_exe_path)
    # Create a list to save the TM-scores
    tmscore_results = []
    # Loop through the fragment_path_list
    for fragment_path in fragment_path_list:
        # Extract the start residue number of the fragment from its filename
        filename = os.path.basename(fragment_path)
        fragment_start = filename[len(query_name):].split("_")[1]
        # Align the fragment with the target structure and save the command standard output
        command_output = tmalign(target_repunit_path, fragment_path)
        # Read and parse the command output to extract the TM-score
        lines = command_output.split("\n")
        tm_score = None
        for line in lines:
            if line.startswith('TM-score=') and 'if normalized by average length' in line:
                # Extract the TM-score value
                tm_score = float(line.split('=')[1].split()[0])
                break  # Stop looping after finding the desired TM-score line
        # Add the start residue number of the fragment (x) and the TM-score (y) as a list tmscore_results
        tmscore_results.append([fragment_start, tm_score])

    # Sort the list based on the start residue numbers
    tm_score_results = [[int(i[0]), i[1]] for i in tmscore_results]
    tm_score_results.sort(key=lambda x: x[0])

    x = [i[0] for i in tm_score_results]
    y = [i[1] for i in tm_score_results]
    # Convert to numpy arrays
    y = np.array(y)
    x = np.array(x)

    return x, y
