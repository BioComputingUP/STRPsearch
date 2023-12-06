# RepeatsDB-lite 2.0
RepeatsDB-lite 2.0 is a specialized tool designed for rapid and precise identification and mapping of structured tandem repeats in proteins (STRPs).

## Getting Started

To get started with the project, you can choose one of the following methods:

### Method 1: Using requirements.txt

1. Install all the dependencies listed in the `requirements.txt` file:
```
   pip install -r requirements.txt
```
Note: Inside the requirements.txt file, you'll find a commented section that includes dependencies which cannot be installed with pip. To install these dependencies, you can use Conda by running the following commands:
```
conda install -c conda-forge -c bioconda foldseek
conda install -c bioconda tmalign
conda install -c conda-forge pymol-open-source
```
2. Run the project with the following command:
```
python3 ./bin/main.py [OPTIONS] COMMAND [ARGS]...
```

### Method 2: Using Conda Environment
1. Import and activate the Conda environment from the `environment.yml` file:
```
conda env create -f environment.yml
conda activate predstrp_env
```
2. Run the project inside the environment with the following command:
```
python3 ./bin/main.py [OPTIONS] COMMAND [ARGS]...
```

### Method 3: Using Docker
1. Build the Docker image using the provided `Dockerfile`:
```
docker build -t repeatsdb-lite .
```
2. To run the container, while being able to interact with the files on the host machine, use the following command:
```
docker run -it --entrypoint /bin/bash -v /home/directory/:/app repeatsdb-lite
```
Be aware that `-v /home/directory/:/app` command mounts the specified directory (`/home/directory/`) to the working directory of the container. This ables the container to read and write files on the host machine.

3. Run the software similar to other methods
```
python3 ./bin/main.py [OPTIONS] COMMAND [ARGS]...
```

## Usage:
The tools has three Commands, each with its positional arguments and options. 

To list the available commands run:

```python3 bin/main.py --help```

Which returns the following commands:

| Command | Description |
|---------|-------------|
| `query-file` | Query an existing PDB/CIF formatted structure file by providing the file path |
| `download-pdb` | Download and query a structure from PDB by providing the PDB ID and the specific Chain of interest |
| `download-model` |  Download and query an AlphaFold model by providing the UniProt ID and the AlphaFold version of interest |
| `version` | Show the version and exit | 

## query-file

### Arguments
* `input_file` (TEXT):  Path to the input structure file to query (PDB/mmCIF). This argument is required. Default: None
* `out_dir` (TEXT): Path to the output directory. This argument is required. Default: None

### Options
* `--chain` (TEXT): Specific chain to query from the structures. Default: all
* `--temp-dir` (TEXT): Path to the temporary directory. Default: /tmp
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4
* `--keep-temp / --no-keep-temp`: Whether to keep the temporary directory and files. Default: no-keep-temp
* `--help`: Show this message and exit

## download-pdb

### Arguments
* `pdb_id` (TEXT): PDB ID of the experimental structure to download and query. This argument is required. Default: None
* `out_dir` (TEXT): Path to the output directory. This argument is required. Default: None

### Options
* `--chain` (TEXT): Specific chain to query from the structures. Default: all
* `--temp-dir` (TEXT): Path to the temporary directory. Default: /tmp
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4
* `--keep-temp / --no-keep-temp`: Whether to keep the temporary directory and files. Default: no-keep-temp
* `--help`: Show this message and exit

## download-model

### Arguments
* `uniprot_id` (TEXT): UniProt ID of the AlphaFold-predicted model to download and query. This argument is required. Default: None
* `af_version` (TEXT): Version of AlphaFold to download predicted models from. This argument is required. Default: None
* `out_dir` (TEXT): Path to the output directory. This argument is required. Default: None

### Options
* `--temp-dir` (TEXT): Path to the temporary directory. Default: /tmp
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4
* `--keep-temp / --no-keep-temp`: Whether to keep the temporary directory and files. Default: no-keep-temp
* `--help`: Show this message and exit

## Examples

If you already have a PDB/CIF formatted structure file and you want to query all the chains in the structure, keeping temporary directory and files:
```
python3 ./bin/main.py query-file /input/file /output/directory --keep-temp
```

If you want to automatically download and query a specific experimental structure from PDB (e.g. chain B of PDB structure 1A0R), without keeping temporary directory and files:
```
python3 ./bin/main.py download-pdb 1a0r /output/directory --chain B
```

If you want to automatically download and query a predicted-model from AlphaFold version 4 (e.g. UniProt ID: Q9HXJ7)
```
python3 ./bin/main.py download-model Q9HXJ7 4 /output/directory 
```


