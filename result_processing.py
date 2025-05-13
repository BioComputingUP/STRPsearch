import os
import json
import csv

def extract_data_from_json(json_path):
    with open(json_path) as f:
        data = json.load(f)

    region_data = next((entry for entry in data if entry["type"] == "region"), None)
    unit_data = [entry for entry in data if entry["type"] == "unit"]
    insertion_data = [entry for entry in data if entry["type"] == "insertion"]

    if not region_data:
        return None

    structure_id = region_data["structure_id"]
    chain_id = region_data["chain_id"]
    cls = region_data["class"]
    topology = region_data["topology"]
    region_start = region_data["start"]
    region_end = region_data["end"]


    units = [f"{unit['start']}-{unit['end']}" for unit in unit_data]
    nb_units = len(units)
    nb_insertions = len(insertion_data)

    return {
        "structure_id": structure_id,
        "chain_id": chain_id,
        "class": cls,
        "topology": topology,
        "region_start": region_start,
        "region_end": region_end,
        "units": ";".join(units),
        "insertions": nb_insertions,
        "nb_units": nb_units
    }

def get_execution_time(root_folder):
    time_file = os.path.join(root_folder, "time.txt")
    if os.path.exists(time_file):
        with open(time_file) as f:
            return f.read().strip()
    return "N/A"

def process_results(root_folder, output_csv):
    exe_time = get_execution_time(root_folder)
    rows = []

    for subdir in os.listdir(root_folder):
        subdir_path = os.path.join(root_folder, subdir)
        if not os.path.isdir(subdir_path):
            continue

        for chain in os.listdir(subdir_path):
            chain_path = os.path.join(subdir_path, chain)
            if not os.path.isdir(chain_path):
                continue

            for region in os.listdir(chain_path):
                region_path = os.path.join(chain_path, region)
                if not os.path.isdir(region_path):
                    continue

                for file in os.listdir(region_path):
                    if file.endswith(".json"):
                        json_path = os.path.join(region_path, file)
                        entry = extract_data_from_json(json_path)
                        if entry:
                            entry["exe_time"] = exe_time
                            rows.append(entry)

    fieldnames = [
        "structure_id", "chain_id", "class", "topology",
        "region_start", "region_end", "units", "insertions", "exe_time", "nb_units"
    ]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[✓] Fichier CSV généré : {output_csv}")


process_results("/home/yusra/Projects/STRPsearch/output/1a0c", "results_summary.csv")
