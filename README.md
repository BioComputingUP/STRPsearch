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
```
2. Run the project with the following command:
```
python3 ./bin/main.py [OPTIONS] COMMAND [ARGS]...
```

### Method 2: Using Conda Environment
1. Import and activate the Conda environment from the `rdblite_env.yml` file:
```
conda env create -f rdblite_env.yml
conda activate rdblite_env
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
2. To run the container, use the following command:
```
docker run -v $(pwd):/app/ repeatsdb-lite analyze-directory input/directory/ output/directory
```
The line `-v $(pwd):/app/` mounts the current directory (`$(pwd)`) to the working directory of the 
container specified in the `Dockerfile` by the line `WORKDIR /app`. This way, the containerized application
is able to access the input and write the output to the local machine.

## Usage:
The tools has three Commands, each with its positional arguments and options. 

To list the available commands run:

```python3 bin/main.py --help```

Which returns the following commands:

| Command | Description |
|---------|-------------|
| `analyze-directory` | Run the pipeline on a directory containing PDB files |
| `download-model` | Run the pipeline by querying a UNIPROT ID and downloading an AlphaFold model |
| `download-pdb` | Run the pipeline downloading a structure and querying a specific chain |
| `version` | Show the version and exit. | 

## Analyze-directory

### Arguments
* `in_dir` (TEXT): Path to directory containing PDB files. This argument is required. Default: None.
* `out_dir` (TEXT): Path to directory where output will be saved. This argument is required. Default: None.

### Options
* `--keep-temp / --no-keep-temp`: Keep temporary files. Default: no-keep-temp.
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01.
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4.
* `--help`: Show this message and exit.

## Download-model

### Arguments
* `uniprot_id` (TEXT): UniProt ID of the AlphaFold structure to query. This argument is required. Default: None.
* `af_version` (TEXT): Version of AlphaFold to download structure from. This argument is required. Default: None.
* `out_dir` (TEXT): Path to directory where output will be saved. This argument is required. Default: None.

### Options
* `--keep-temp` / `--no-keep-temp`: Keep temporary files. Default: no-keep-temp.
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01.
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4.
* `--help`: Show this message and exit.

## Download-pdb

### Arguments
* `pdb_id` (TEXT): PDB ID to download. This argument is required. Default: None.
* `pdb_chain` (TEXT): PDB chain to query. This argument is required. Default: None.
* `out_dir` (TEXT): Path to directory where output will be saved. This argument is required. Default: None.

### Options
* `--keep-temp` / `--no-keep-temp`: Keep temporary files. Default: no-keep-temp.
* `--max-eval` (FLOAT): Maximum E-value of the targets to prefilter. Default: 0.01.
* `--min-height` (FLOAT): Minimum height of TM-score signals to be processed. Default: 0.4.
* `--help`: Show this message and exit.

## Examples

If you already have one or more pdb/mmcif format structures of **single polypeptide chains** stored in a directory, while also keeping temporary files:
```
python3 ./bin/main.py analyze-directory /input/directory -o /output/directory --keep-temp 
```

If you want to download a specific structure from PDB (e.g. chain C of 4g8l PDB structure), without keeping temporary files:
```
python3 ./bin/main.py download-pdb 4g8l C /output/directory 
```

If you want to download a predicted structure from AlphaFold (e.g. UniProt ID: Q05823)
```
python3 ./bin/main.py download-model Q05823 4 /output/directory 
```


