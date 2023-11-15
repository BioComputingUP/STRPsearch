from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from datetime import datetime
from Bio.PDB import Select
import numpy as np
import json


def get_res_index(res_num: str, chain_residues: list):
    """
    Returns the residue index of a residue
    """

    # Loops through the residues of a given list of residues
    # Finds the index of the specified residue number in the list
    # (Not always the index and the residue number are the same e.g. chimeric proteins)
    res_index = [i for i, chain_res in enumerate(chain_residues) if str(chain_res.id[1]) == str(res_num)]
    return res_index[0]


def get_chain_range(start, end, chain_residues: list):
    """
    Returns the range of subsection in a chain
    """

    # Get the residue indices for the start and end residues on the given list of residues
    chain_start = get_res_index(start, chain_residues)
    chain_end = get_res_index(end, chain_residues)
    # Slice and return that range of the residue list (end residue inclusive)
    chain_range = chain_residues[chain_start:chain_end + 1]
    return chain_range


def get_res_frames(res_range, length, step):
    """
    Returns a list of fragmented frames of a residue range
    """

    # Create a list to save the fragment frames
    frames = []
    # Define the start index
    start = 0
    # Loop for length of the range
    for i in range(len(res_range)):
        # Slice the range by indices and save it as frame
        # The end of each frame is defined as the start + the defined length of the frame
        frame = res_range[start:start + length]
        # Do it until the length of the frames fall no shorter that 75% of the defined length
        if len(frame) >= round(length * 0.75):
            # Add the frame to the list of fragment frames
            frames.append(frame)
            # Increment the start by the defined step
            start += step
    return frames


def get_structure(res_range, chain_letter, structure, out_path, io_handler):
    """T
    rims and saves the structure in the range of start and end residues on a chain
    """

    # Select residues within the specified range and chain
    class ResSelect(Select):
        def accept_residue(self, res):
            if res in res_range and res.parent.id == chain_letter:
                return True
            else:
                return False

    # Create a PDBIO object, set the structure, and save the selected residues to the output path
    io_handler.set_structure(structure)
    io_handler.save(out_path, ResSelect())


def calculate_ranges(peak_residue_nums, distance, max_length, flexibility):
    """
    Defines the ranges of repeat units by the given residues that corresponds to TM-score graph peaks
    """

    # Define the maximum allowed distance between the units
    # Simply the expected distance + a certain percentage of flexibility in length
    max_distance = round(distance + distance * flexibility)
    # Create a list to save the ranges
    ranges = []
    # Loop for the length of the given list of peak residues
    for i in range(len(peak_residue_nums)):
        # Get the res number corresponding to the index
        current_num = peak_residue_nums[i]
        # Check if the current peak number is not the last number in the list
        if i + 1 < len(peak_residue_nums):
            # If not, get the next number in the list
            next_num = peak_residue_nums[i + 1]
            # Check if the distance to the next peak is less than or equal to the max allowed distance
            if next_num - current_num <= max_distance:
                # If yes, the end of the range would be the next peak position
                start_range = current_num
                end_range = next_num - 1
                ranges.append((start_range, end_range))
            else:
                # If not, the end of the range would be the start + the expected distance
                start_range = current_num
                end_range = current_num + distance
                ranges.append((start_range, end_range))
        else:
            # If the current peak res number is the last one in the list
            start_range = current_num
            # Check if plus the usual distance, it falls shorter than the max length of the chain
            if (current_num + distance) <= max_length:
                # If yes, define the range accordingly
                end_range = current_num + distance
            # If not, involve the allowed flexibility and see if it falls shorter than the max length of the chain
            elif current_num + round(distance * (1 - flexibility)) <= max_length:
                end_range = max_length
            # If not any cases above, ignore and break the loop
            else:
                break
            # Add the defined range to the list
            ranges.append((start_range, end_range))
    return ranges


def make_json(pdb_id_chain, ct, regions_dict, out_path):
    """
    Creates and saves a JSON file containing the annotation of the repeat region/regions of a specific type
    """

    # Create a list to save the annotation dictionaries
    entries_list = []

    # Extract the essential data from the given arguments
    pdb_id = pdb_id_chain[:-1]
    pdb_chain = pdb_id_chain[-1]
    class_type = ct.split(".")[0]
    topology_type = ct.split(".")[1]

    # Loop through the regions and their annotation in the regions_dict
    # Embed the essential data with the specified format
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

    # Dump the list of dicts into a JSON
    json_dict = json.dumps(entries_list)
    # Save at the specified path
    with open(out_path, 'w') as fp:
        fp.write(json_dict)


def parse_regions(pdb_id_chain, units_list, insert_len):
    """
    Based on the number, position, and the range of the predicted units, parse one or more repeat regions
    """

    # Create a dict to save the units associated with each region indicated by numbers
    primary_regions_dict = {1: []}
    # Variable counts associates with the number of mapped regions
    count = 1

    # Loop through the units_list
    for index, unit in enumerate(units_list):
        # Check if the unit is not the last unit
        if index + 1 < len(units_list):
            # If not, calculate the distance between the end residue of the unit and start residue of the next unit
            inter_unit_distance = units_list[index + 1][0] - unit[1]
            # Check if it is smaller than or equal to the allowed length of insertion
            if inter_unit_distance <= insert_len:
                # If it is, append the unit to the list associated with the region count in the primary_regions_dict
                primary_regions_dict[count].append(unit)
            else:
                # If not, append the unit to the list with the current region count
                primary_regions_dict[count].append(unit)
                # Increment the count by one, so it maps an extra region from this point
                count += 1
                primary_regions_dict[count] = []
        else:
            # If it is the last unit, just added to the list with the current region count
            primary_regions_dict[count].append(unit)

    # Create a dict to save the formatted repeat regions
    final_regions_dict = {}
    # Loop through the primary_regions_dict
    for region_num, ranges in primary_regions_dict.items():
        # If the region contains more than potential units
        if len(ranges) >= 3:
            # Define the region_id similar to PDB-id-chain_region-start-res_region-end-res
            region_id = f"{pdb_id_chain}_{ranges[0][0]}_{ranges[-1][1]}"
            # Exactly define the range of units and insertions and save into their lists accordingly
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
            # Add the integral components as a dict to the final dict with the key of their region_id
            final_regions_dict[region_id] = {"units": unit_ranges, "insertions": insertion_ranges}

    return final_regions_dict


def adjust_graph_ends(x, y, frame_step=1):
    """
    Slightly lowers the y of the first and last point in the graph
    """

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
    """
    Plots and saves a TM-score graph highlighted by integral components of the repeat region/regions
    """

    # Define the colors for alternate coloring of units
    colors = ['blue', 'red']
    fig, ax = plt.subplots(dpi=200)

    # Plot the line graph
    ax.plot(x, y)
    # Plot the dots associated with residue numbers on the graph
    ax.scatter(x, y, marker='o', color='blue', s=2)
    # Loop their repeat regions in regions_dict and highlight the components on the plot
    for region_id, v in regions_dict.items():
        for i, (start, end) in enumerate(v["units"]):
            color = colors[i % len(colors)]
            ax.axvspan(start, end, facecolor=color, alpha=0.3)

        for i, (start, end) in enumerate(v["insertions"]):
            ax.axvspan(start, end, facecolor="yellow", alpha=0.3)

    # Format the plot
    ax.set_xlabel("Residue Number")
    ax.set_ylabel("TM-score")
    plt.tight_layout()
    # Save at the specified path in png format
    plt.savefig(out_path, format="png")


def smooth_graph(y, target_avg_len, window_p):
    """
    Smooths a graph data (heights)
    """

    # Define the window size based the provided parameters
    window_size = round(target_avg_len * window_p)
    if window_size % 2 == 0:
        window_size -= 1

    poly_order = 0
    # Smooth the graph by savgol_filter module
    y = savgol_filter(y, window_size, poly_order)

    return y

def get_current_time():
    """
    Returns the current time in the format of HH:MM:SS
    """

    now = datetime.now()
    formatted_time = now.strftime("%H:%M:%S")
    return formatted_time
