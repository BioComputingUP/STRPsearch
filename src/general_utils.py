from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from datetime import datetime
from Bio.PDB import Select
from Bio import BiopythonWarning
import seaborn as sns
import numpy as np
import json, os, math, warnings, gemmi, re
from Bio.PDB import PDBParser, MMCIFIO , PDBIO,MMCIFParser
from contextlib import redirect_stdout

warnings.filterwarnings("ignore", category=BiopythonWarning)


from protein_domain_segmentation import ChainsawCluster
class ResidueRangeSelect(Select):
    def __init__(self, chain_id, start_res, end_res):
        self.chain_id = chain_id
        self.start_res = start_res
        self.end_res = end_res

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (
            residue.get_parent().id == self.chain_id and
            self.start_res <= residue.id[1] <= self.end_res
        )

def extract_segment_to_cif(pdb_file, chain_id, start_res, end_res, output_file):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        # Check chain exists
        chains = [c.id for c in structure[0]]
        if chain_id not in chains:
            print(f"Warning: chain {chain_id} not found in {pdb_file}")
            return False
        
        # Check residues exist
        residues = [r for r in structure[0][chain_id] if start_res <= r.id[1] <= end_res]
        if not residues:
            print(f"Warning: no residues in range {start_res}-{end_res} for chain {chain_id} in {pdb_file}")
            return False
        
        # Save CIF
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(output_file, select=ResidueRangeSelect(chain_id, start_res, end_res))
        return True
    
    except Exception as e:
        print(f"Error saving segment {start_res}-{end_res} for chain {chain_id} in {pdb_file}: {e}")
        return False


def get_chain_id_from_filename(filename):
    """
    Extracts the chain ID from a .cif filename in the format '*_chain.cif',
    taking the last underscore-separated part as the chain ID.

    Args:
        filename (str): The name of the .cif file (e.g., '1a0t_P.cif' or 'some_file_name_P.cif').

    Returns:
        str: The chain ID extracted from the filename.
    """
    base_name = os.path.splitext(filename)[0]  # Remove extension
    if "_" not in base_name:
        raise ValueError(f"Filename '{filename}' does not contain a chain ID.")
    chain_id = base_name.split("_")[-1]  # Take last part
    return chain_id



def extract_regions(region_string):
    """
    Extracts start and end values from a region string and returns a list of dictionaries.
    Handles cases where residue numbers may have insertion codes (e.g., 365I).
    """
    if not region_string:
        return [] 

    regions = region_string.split(",")
    result = []
    for region in regions:
        regions_l = region.split("_")
        for reg in regions_l:
            parts = [p for p in reg.split("-") if p]  # Remove empty strings
            if len(parts) >= 2:
                # Extract only the numeric part for start and end
                def extract_num(s):
                    m = re.match(r"(\d+)", s)
                    return int(m.group(1)) if m else None
                start = extract_num(parts[0])
                end = extract_num(parts[1])
                if start is not None and end is not None:
                    result.append({"start": start, "end": end})
                else:
                    print(f"Skipping malformed region: {reg}")
            else:
                print(f"Skipping malformed region: {reg}")

    return result

def is_polymer_chain_cif(filepath):
    """
    Détermine si un fichier .cif correspond à une chaîne protéique (ou un polymère biologique).
    Retourne True s'il s'agit d'une chaîne (avec des résidus comme ALA, GLY...), False sinon.
    """
    try:
        with open(filepath, 'r') as file:
            has_residues = False
            for line in file:
                # Cherche les lignes ATOM (et pas HETATM)
                if line.startswith("ATOM"):
                    # Vérifie si le résidu appartient à une chaîne classique d'acides aminés
                    tokens = line.split()
                    if len(tokens) > 5:
                        residue = tokens[5]
                        if residue in {
                            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                            "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                            "PRO", "SER", "THR", "TRP", "TYR", "VAL"
                        }:
                            has_residues = True
                            break
            return has_residues
    except Exception as e:
        print(f"Error : {e}")
        return False



def segment_cif_directory(input_dir, output_dir):
    """
    Processes a directory of .cif files, applies Chainsaw to predict chopping regions,
    extracts segments, and deletes temporary .pdb files.

    Args:
        input_dir (str): Path to the directory containing .cif files.
        output_dir (str): Path to the directory where output .cif files will be saved.
    """
    # os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        chainsaw_cluster = ChainsawCluster()  # Initialize ChainsawCluster
        for cif_file in os.listdir(input_dir):
            
            if cif_file.endswith(".cif"):
                cif_path = os.path.join(input_dir, cif_file)
                if is_polymer_chain_cif(cif_path):
                    chain_id = get_chain_id_from_filename(cif_file)  # Extract chain ID from filename
                    pdb_file = cif_path.replace(".cif", ".pdb")  # Temporary .pdb file path

                    # Convert .cif to .pdb
                    parser = MMCIFParser(QUIET=True)
                    
                    try:
                        structure = parser.get_structure("cif_file", cif_path)
                    except IndexError as e:
                        print(f"IndexError while parsing {cif_file}: {e}")
                        return False  # or handle as appropriate
                    except KeyError as e:
                        print(f"KeyError while parsing {cif_file}: {e}")
                        return False
                    except Exception as e:
                        print(f"Unexpected error while parsing {cif_file}: {e}")
                        return False
                    for model in structure:
                        chains_to_remove = [chain for chain in model if len(chain.id) > 1]
                        for chain in chains_to_remove:
                            print(f"Skipping chain id '{chain.id}' (invalid for PDB format).")
                            model.detach_child(chain.id)
                    io = PDBIO()
                    io.set_structure(structure)
                    io.save(pdb_file)
                    # Apply Chainsaw to predict chopping regions
                    with open(os.devnull, 'w') as devnull:
                        with redirect_stdout(devnull):
                            try:
                                chainsaw_result = chainsaw_cluster.predict_from_pdb(pdb_file)
                            except RuntimeError as e:
                                print(f"RuntimeError during Chainsaw prediction for {pdb_file}: {e}")
                                continue  # Skip this file and continue with the next
                            except Exception as e:
                                print(f"Unexpected error during Chainsaw prediction for {pdb_file}: {e}")
                                continue
                    # Extract regions from Chainsaw results
                    regions = extract_regions(chainsaw_result)

                    # Extract segments for each region
                    if not regions:
                        print(f"no segmentation found for {cif_file}")
                    else:
                        for region in regions:
                            start = region['start']
                            end = region['end']
                            output_cif = os.path.join(
                                output_dir, f"{start}_{end}_{os.path.splitext(cif_file)[0]}.cif"
                            )
                            extract_segment_to_cif(pdb_file, chain_id, start, end, output_cif)
                        os.remove(cif_path)  # Remove the original .cif file

                    # Delete the temporary .pdb file
                    os.remove(pdb_file)


def get_res_index(res_num: str, chain_residues: list):
    """
    Returns the residue index of a residue
    """

    # Loops through the residues of a given list of residues
    # Finds the index of the specified residue number in the list
    # (Not always the index and the residue number are the same e.g. chimeric proteins)
    res_index = [i for i, chain_res in enumerate(chain_residues) if str(chain_res.seqid.num) == str(res_num)]
    return res_index[0]


def find_largest_smaller_number(res_num: int, chain_residues: list):
    """
    Finds the largest residue number in the chain that is smaller than or equal to the given residue number.
    """
    smaller_numbers = [int(chain_res.seqid.num) for chain_res in chain_residues if int(chain_res.seqid.num) <= res_num]
    largest_smaller_number = max(smaller_numbers)
    return largest_smaller_number


def get_chain_range(start, end, chain_residues: list):
    """
    Returns a list of residue numbers (integers) in the specified range of a gemmi chain.

    Args:
        start (int): The starting residue number.
        end (int): The ending residue number.
        chain_residues (list): A list of gemmi residues.

    Returns:
        list: A list of residue numbers (integers) in the specified range.
    """
    # Get the residue indices for the start and end residues
    chain_start = get_res_index(start, chain_residues)
    chain_end = get_res_index(end, chain_residues)

    # Slice the range of residues and extract their residue numbers
    chain_range = [res.seqid.num for res in chain_residues[chain_start:chain_end + 1]]
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
        # Do it until the length of the frames fall no shorter than 75% of the defined length
        if len(frame) >= round(length * 0.75):
            # Add the frame to the list of fragment frames
            frames.append(frame)
            # Increment the start by the defined step
            start += step
    return frames


def get_structure(res_range, chain_letter, structure, out_path, io_handler):
    """
    Trims and saves the structure in the range of start and end residues on a chain
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


def get_structure_cif(res_range, chain_id, structure: gemmi.Structure, out_path: str):
    """
    Trims and saves a .cif structure using gemmi.
    res_range: list of integers (residue numbers)
    """
    model = structure[0]
    original_chain = model[chain_id]

    # Create a new structure with the selected residues
    new_structure = gemmi.Structure()
    new_structure.name = structure.name
    new_model = gemmi.Model('1')
    new_chain = gemmi.Chain(chain_id)

    for res in original_chain:
        if res.seqid.num in res_range:
            new_chain.add_residue(res.clone())

    new_model.add_chain(new_chain)
    new_structure.add_model(new_model)

    # Convert the new structure into a CIF block
    cif_block = new_structure.make_mmcif_block()

    # Convert the CIF block to a string
    cif_string = cif_block.as_string()

    # Write the CIF string to the output file
    with open(out_path, 'w') as cif_file:
        cif_file.write(cif_string)


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


def make_json(structure_id, chain_id, ct, region_id, components, out_path):
    """
    Creates and saves a JSON file containing the annotation of the repeat region/regions of a specific type
    """

    # Create a list to save the annotation dictionaries
    entries_list = []

    # Extract the essential data from the given arguments
    class_type = ct.split(".")[0]
    topology_type = ct.split(".")[1]

    # Loop through the regions and their annotation in the regions_dict
    # Embed the essential data with the specified format
    for unit in components["units"]:
        start = str(unit[0])
        end = str(unit[1])
        unit_dict = {
            "structure_id": structure_id,
            "chain_id": chain_id,
            "type": "unit",
            "start": start,
            "end": end,
            "class": class_type,
            "topology": topology_type,
            "origin": "Predicted",
            "reviewed": False,
            "region_id": region_id
        }

        entries_list.append(unit_dict)

    for insertion in components["insertions"]:
        start = str(insertion[0])
        end = str(insertion[1])
        insertion_dict = {
            "structure_id": structure_id,
            "chain_id": chain_id,
            "type": "insertion",
            "start": start,
            "end": end,
            "class": class_type,
            "topology": topology_type,
            "origin": "Predicted",
            "reviewed": False,
            "region_id": region_id
        }

        entries_list.append(insertion_dict)

    region_dict = {
        "structure_id": structure_id,
        "chain_id": chain_id,
        "type": "region",
        "start": f"{components['units'][0][0]}",
        "end": f"{components['units'][-1][1]}",
        "class": class_type,
        "topology": topology_type,
        "origin": "Predicted",
        "reviewed": False,
        "region_id": region_id
    }

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
    # Compute the new value as less than the current last value
    y_new_value2 = round(y_last_value * 0.99, 5)
    x_new_value2 = round(x_last_value + frame_step)
    # Create a new array with the new values and concatenate it with the existing array
    y = np.concatenate(([y_new_value], y))
    y = np.concatenate((y, [y_new_value2]))
    x = np.concatenate(([x_new_value], x))
    x = np.concatenate((x, [x_new_value2]))
    return x, y


def plot_tmscore_graph(x, y, region_components, out_path,json_out_path):
    from matplotlib.lines import Line2D

    """
    Plots and saves a TM-score graph highlighted by integral components of the repeat region/regions
    """

    # Set a Seaborn style
    sns.set_style("ticks")
    sns.despine()
    # Define the colors for alternate coloring of units
    colors = ["salmon", "cornflowerblue"]

    # Plot the line graph
    fig, ax = plt.subplots(dpi=200)
    ax.plot(x, y, color="black", linewidth=1.2)

    # Loop their repeat regions in regions_dict and highlight the components on the plot
    for i, (start, end) in enumerate(region_components["units"]):
        color = colors[i % len(colors)]
        ax.axvspan(start, end, facecolor=color, alpha=0.85)

    for i, (start, end) in enumerate(region_components["insertions"]):
        ax.axvspan(start, end, facecolor="khaki", alpha=0.85)

    # Format the plot
    ax.set_xlabel("Residue Number")
    ax.set_ylabel("TM-score")

    # plt.xticks(fontsize=16)
    # plt.yticks(fontsize=16)
    plt.tight_layout()
    # Save at the specified path in png format
    plt.savefig(out_path, format="png")
    plt.close(fig)
    # clean_units = [
    #     [int(start), int(end)] for start, end in region_components["units"]
    # ]
    # clean_insertions = [
    #     [int(start), int(end)] for start, end in region_components["insertions"]
    # ]

    # json_data = {
    #     "x": x.tolist() if isinstance(x, np.ndarray) else list(x),
    #     "y": y.tolist() if isinstance(y, np.ndarray) else list(y),
    #     "regions": {
    #         "units": clean_units,
    #         "insertions": clean_insertions 
    #     }
    # }
    
    # with open(json_out_path, "w") as f:
    #     json.dump(json_data, f, indent=4)



def smooth_graph(y, target_avg_len, window_p):
    """
    Smooths a graph data (heights)
    """

    # Define the window size based the provided parameters
    window_size = round(target_avg_len * window_p)
    if window_size % 2 == 0:
        window_size -= 1
    # Ensure window_size is at least 3 and not greater than len(y)
    window_size = max(3, min(window_size, len(y) if len(y) % 2 == 1 else len(y) - 1))
    # Ensure window_size is at least 3 and not greater than len(y)
    window_size = max(3, min(window_size, len(y) if len(y) % 2 == 1 else len(y) - 1))

    poly_order = 0
    # Smooth the graph by savgol_filter module
    y = savgol_filter(y, window_size, poly_order)

    return y



def export_graph_data(x, y, region_components, out_path):
    """
    Exports plot data to JSON for the Angular interactive UI
    """
    data = {
        "residues": x.tolist() if hasattr(x, 'tolist') else list(x),
        "tm_scores": y.tolist() if hasattr(y, 'tolist') else list(y),
        "units": region_components["units"],       # List of [start, end]
        "insertions": region_components["insertions"] # List of [start, end]
    }
    
    with open(out_path.replace(".png", ".json"), "w") as f:
        json.dump(data, f)

def time():
    """
    Returns the current time in the format of HH:MM:SS
    """

    now = datetime.now()
    formatted_time = now.strftime("%H:%M:%S")
    return formatted_time


def create_pymol_session(region_id, structure_path, components, output_path):
    """
    Creates a PyMOL session that the repeat components are highlighted on the associated structure
    """
    from pymol import cmd as pcmd
    # Load the structure
    pcmd.load(structure_path, region_id)

    # Loop the units from the components and highlight them alternatively
    for idx, (start, end) in enumerate(components["units"]):
        selection_name = f'unit_{idx + 1}'
        pcmd.select(selection_name, f'resi {start}-{end}')
        color_name = 'raspberry' if idx % 2 == 0 else 'marine'  # Alternating colors
        pcmd.color(color_name, selection_name)

    # Loop the insertions from components and highlight them in yellow
    for idx, (start, end) in enumerate(components["insertions"]):
        selection_name = f"insertion_{idx + 1}"
        pcmd.select(selection_name, f'resi {start}-{end}')
        color_name = 'paleyellow'
        pcmd.color(color_name, selection_name)

    # Save the PyMOL session
    pcmd.save(output_path)
    # Reinitialize the PyMOL for future runs
    pcmd.reinitialize()


def check_files(in_dir):
    """
    Checks filenames with pdb/cif extension in a directory,
    returns two lists:
    1. All the filepaths with pdb/cif extension
    2. The ones with spacing in their filenames
    """

    # Create the two lists
    pdb_cif_fps = []
    have_space_fps = []
    # Loop through the filenames of the directory
    for filename in os.listdir(in_dir):
        # Check if the filename is PDB or mmCIF
        if filename.endswith(".pdb") or filename.endswith(".cif"):
            # Define the filepath
            file_path = os.path.join(in_dir, filename)
            # Add to the general list
            pdb_cif_fps.append(file_path)
            # If space in the filename, add to the spacing list
            if " " in filename:
                have_space_fps.append(file_path)
    return pdb_cif_fps, have_space_fps


def get_repeat_classi(ontology_df, code):
    """
    Retrieve the classification name of the repeat with its Class. Topology code
    """
    if code is None or (isinstance(code, float) and math.isnan(code)):
        # Optionally: return None or a default string instead of raising
        return None
    filtered = ontology_df[ontology_df["Code"] == code]
    if filtered.empty:
        return None
    return filtered["Name"].values[0]
