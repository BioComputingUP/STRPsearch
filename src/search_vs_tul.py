from src.foldseek import Foldseek
import pandas as pd


def search_tul(foldseek_exe_path, query_dir, tul_fs_db, output_file, temp_dir):
    foldseek = Foldseek(exe_path=foldseek_exe_path)
    columns = "query,target,evalue,qstart,qend"
    foldseek.easy_search(input_file=query_dir, db=tul_fs_db, output_file=output_file, temp_dir=temp_dir,
                         columns=columns)


def check_overlap(range1, range2):
    range1 = set(range1)
    overlapped = len(range1.intersection(range2))
    overlap_percentage = overlapped / len(range1)
    if overlap_percentage >= 0.5:
        return True
    return False


def filter_overlap(df):
    master_index = []

    min_eval_index = df["e_value"].idxmin()
    master_index.append(min_eval_index)
    for q_i in range(len(df)):
        if q_i not in master_index:
            q_row = df.iloc[q_i]
            q_eval = q_row["e_value"]
            q_range = range(int(q_row["q_start"]), int(q_row["q_end"]))
            overlap = []
            for idx, master_i in enumerate(master_index):
                master_row = df.iloc[master_i]
                master_eval = master_row["e_value"]
                master_range = range(int(master_row["q_start"]), int(master_row["q_end"]))

                if check_overlap(master_range, q_range):
                    overlap.append(True)
                    if q_eval < master_eval:
                        master_index[idx] = q_i
                else:
                    overlap.append(False)

            if True not in overlap:
                master_index.append(q_i)

    return df.iloc[master_index]


def find_target(output_file, max_eval):
    df = pd.read_csv(output_file, delimiter="\t", header=None, dtype={"ctr": str})
    df.columns = ["query", "target", "e_value", "q_start", "q_end"]
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

    target_list = []
    for query in df["query"].unique():
        query_df = df[df["query"] == query]
        # Group the DataFrame by unique class topology values and find the index with the lowest eval for each group
        lowest_eval_indices = query_df.groupby('t_ct')["e_value"].idxmin()

        # Extract the rows with the lowest eval for each class topology
        lowest_eval_rows = query_df.loc[lowest_eval_indices]
        for i in lowest_eval_rows.index:
            row = lowest_eval_rows.loc[i]
            target_list.append({"query": row["query"], "target": row["target"], "q_start": row["q_start"],
                                "q_end": row["q_end"], "t_ct": row["t_ct"], "t_avg_length": row["t_avg_length"],
                                "e_value": row["e_value"]})

    target_df = pd.DataFrame(target_list)
    target_df = filter_overlap(target_df)

    return True, target_df

