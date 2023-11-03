import json
from Bio.PDB import Select
import numpy as np
import matplotlib.pyplot as plt


def get_res_index(res_num: str, chain_residues: list):
    """Returns the residue index of a residue"""
    res_index = [i for i, chain_res in enumerate(chain_residues) if str(chain_res.id[1]) == str(res_num)]
    return res_index[0]


def get_chain_range(start, end, chain_residues: list):
    """Returns the range of subsection in a chain"""
    chain_start = get_res_index(start, chain_residues)
    chain_end = get_res_index(end, chain_residues)
    chain_range = chain_residues[chain_start:chain_end + 1]
    return chain_range


def get_res_frames(res_range, length, step):
    """Returns a list of fragmented frames of a residue range"""
    frames = []
    start = 0
    for i in range(len(res_range)):
        frame = res_range[start:start + length]
        if len(frame) >= round(length * 0.75):
            frames.append(frame)
            start += step
    return frames


def get_structure(res_range, res_chain, structure, out_path, io_handler):
    """Trims and saves the structure in the range of start and end residues on a chain"""

    # Select residues within the specified range and chain
    class ResSelect(Select):
        def accept_residue(self, res):
            if res in res_range and res.parent.id == res_chain:
                return True
            else:
                return False

    # Create a PDBIO object, set the structure, and save the selected residues to the output path
    io_handler.set_structure(structure)
    io_handler.save(out_path, ResSelect())


def calculate_ranges(numbers, distance, max_length, flexibility):
    """Defines unit ranges"""
    ranges = []
    for i in range(len(numbers)):
        current_num = numbers[i]
        if i + 1 < len(numbers):
            next_num = numbers[i + 1]
            max_distance = round(distance + (distance * flexibility))
            if next_num - current_num <= max_distance:
                start_range = current_num
                end_range = next_num - 1
                ranges.append((start_range, end_range))
            else:
                start_range = current_num
                end_range = current_num + distance
                ranges.append((start_range, end_range))
        else:
            start_range = current_num
            if (current_num + distance) <= max_length:
                end_range = current_num + distance
            elif round(current_num + distance * flexibility) <= max_length:
                end_range = max_length
            else:
                break
            ranges.append((start_range, end_range))
    return ranges


def make_json(pdb_id_chain, ct, regions_dict, out_path):
    entries_list = []
    pdb_id = pdb_id_chain[:-1]
    pdb_chain = pdb_id_chain[-1]
    class_type = ct.split(".")[0]
    topology_type = ct.split(".")[1]
    for region_id, v in regions_dict.items():
        for unit in v["units"]:
            start = str(unit[0])
            end = str(unit[1])
            unit_dict = {"start": start, "end": end, "type": "unit", "pdb_id": pdb_id, "pdb_chain": pdb_chain,
                         "repeatsdb_id": pdb_id_chain, "class": class_type, "topology": topology_type, "fold": "0",
                         "clan": "0", "class_topology": ct, "class_topology_fold": f"{ct}.0",
                         "class_topology_fold_clan": f"{ct}.0.0", "origin": "Predicted", "reviewed": False,
                         "annotator": "", "region_id": region_id, "region_units_num": 0,
                         "region_average_unit_length": 0}

            entries_list.append(unit_dict)

        for insertion in v["insertions"]:
            start = str(insertion[0])
            end = str(insertion[1])
            insertion_dict = {"start": start, "end": end, "type": "insertion", "pdb_id": pdb_id, "pdb_chain": pdb_chain,
                              "repeatsdb_id": pdb_id_chain, "class": class_type, "topology": topology_type, "fold": "0",
                              "clan": "0", "class_topology": ct, "class_topology_fold": f"{ct}.0",
                              "class_topology_fold_clan": f"{ct}.0.0", "origin": "Predicted", "reviewed": False,
                              "annotator": "", "region_id": region_id, "region_units_num": 0,
                              "region_average_unit_length": 0}
            entries_list.append(insertion_dict)

        region_dict = {"start": f"{v['units'][0][0]}", "end": f"{v['units'][-1][1]}", "type": "region",
                       "pdb_id": pdb_id, "pdb_chain": pdb_chain, "repeatsdb_id": pdb_id_chain, "class": class_type,
                       "topology": topology_type, "fold": "0", "clan": "0", "class_topology": ct,
                       "class_topology_fold": f"{ct}.0", "class_topology_fold_clan": f"{ct}.0.0",
                       "origin": "Predicted", "reviewed": False, "annotator": "", "region_id": region_id,
                       "region_units_num": 0, "region_average_unit_length": 0}

        entries_list.append(region_dict)

    json_dict = json.dumps(entries_list)
    with open(out_path, 'w') as fp:
        fp.write(json_dict)


def parse_regions(pdb_id_chain, units_list, insert_len):
    a_dict = {1: []}
    count = 1

    for index, unit in enumerate(units_list):
        if index + 1 < len(units_list):
            if units_list[index + 1][0] - unit[1] <= insert_len:
                a_dict[count].append(unit)
            else:
                a_dict[count].append(unit)
                count += 1
                a_dict[count] = []
        else:
            a_dict[count].append(unit)

    b_dict = {}

    for key, ranges in a_dict.items():
        if len(ranges) >= 3:
            new_key = f"{pdb_id_chain}_{ranges[0][0]}_{ranges[-1][1]}"
            unit_ranges = []
            insertion_ranges = []
            for i in range(len(ranges) - 1):
                current_range = ranges[i]
                next_range = ranges[i + 1]
                unit_ranges.append(current_range)
                diff_start = current_range[1] + 1
                diff_end = next_range[0] - 1
                if diff_start <= diff_end and diff_end - diff_start >= 3:
                    insertion_ranges.append((diff_start, diff_end))
            unit_ranges.append(ranges[-1])
            b_dict[new_key] = {"units": unit_ranges, "insertions": insertion_ranges}

    return b_dict


def adjust_graph_ends(x, y, frame_step=1):
    """Slightly lowers the y of the first and last point in the graph"""
    # Extract the first value from the existing array
    y_first_value = y[0]
    x_first_value = x[0]
    # Extract the last value from the existing array
    y_last_value = y[-1]
    x_last_value = x[-1]
    # Compute the new value as less than the current first value
    y_new_value = round(y_first_value * 0.99, 5)
    x_new_value = round(x_first_value - frame_step)
    # Compute the new value as less than the current first value
    y_new_value2 = round(y_last_value * 0.99, 5)
    x_new_value2 = round(x_last_value + frame_step)
    # Create a new array with the new value and concatenate it with the existing array
    y = np.concatenate(([y_new_value], y))
    y = np.concatenate((y, [y_new_value2]))
    x = np.concatenate(([x_new_value], x))
    x = np.concatenate((x, [x_new_value2]))
    return x, y


def plot_tmscore_graph(x, y, regions_dict, out_path):
    colors = ['red', 'blue']
    fig, ax = plt.subplots(dpi=200)

    ax.plot(x, y)
    ax.scatter(x, y, marker='o', color='blue', s=2)

    for region_id, v in regions_dict.items():
        for i, (start, end) in enumerate(v["units"]):
            color = colors[i % len(colors)]
            ax.axvspan(start, end, facecolor=color, alpha=0.3)

        for i, (start, end) in enumerate(v["insertions"]):
            ax.axvspan(start, end, facecolor="yellow", alpha=0.3)

    ax.set_xlabel("Residue Number")
    ax.set_ylabel("TM-score")
    plt.tight_layout()
    plt.savefig(out_path, format="png")
