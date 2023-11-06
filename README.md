# Bugfix/Refactoring branch

# RepeatsDB-lite 2.0
RepeatsDB-lite 2.0 is a specialized tool designed for rapid and precise identification and mapping of structured tandem repeats in proteins (STRPs).

## Bugfix/Checklist

- [ ] Fix usage
- [ ] Fix installation with pip
- [x] Include -h and -v options
- [x] Remove need to specify if structure exists
- [ ] Make keep tempfiles optional
- [x] Implement functionality with supplied directory
- [x] Implement functionality with PDB download
- [ ] Implement functionality with AlphaFold download

## New bugs found
- Supplying .cif files returns ```The query file format is ambiguous for query XXXX.cif_U```. The added **_U** I do not know why or how it happens.

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
python3 ./bin/repeatsdb-lite-2.py [OPTIONS]
```

### Method 2: Using Conda Environment
1. Import and activate the Conda environment from the `rdblite_env.yml` file:
```
conda env create -f rdblite_env.yml
conda activate rdblite_env
```
2. Run the project inside the environment with the following command:
```
python3 ./bin/repeatsdb-lite-2.py [OPTIONS]
```

### Method 3: Using Docker
1. Build the Docker image using the provided `Dockerfile`:
```
docker build -t repeatsdb-lite .
```
2. Run the Docker container with the following command:
```
docker run repeatsdb-lite [OPTIONS]
```
Note: Please be aware that the Docker container will not have access to paths located on the host system. To provide input and retrieve output, you'll need to transfer files between the container and the host manually.

## List of Options:
| Option | Type	| Description |
| --- | --- | --- |
| --structure_exists (-e) | str | Structure already exists at input directory (1) or not (0) |
| --out_dir (-o) | str | Path to the output directory |
| --temp_dir (-t) | str | Path to the temporary directory |
| --in_dir (-i) | str (Default: None) | Path to the input directory |
| --pdb_id | str (Default: None) | PDB structure to download |
| --pdb_chain | str (Default: None) | PDB chain to query |
| --uniprot_id | str (Default: None) | UniProt ID of the AlphaFold structure to query |
| --af_version | str (Default: None) | Version of AlphaFold to download structure from |
| --max_eval | float (Default: 0.01) | Maximum E-value of the targets to prefilter |
| --min_height | float (Default: 0.4) | Minimum height of TM-score signals to be processed |


## Examples

If you already have one or more pdb/mmcif format structures of **single polypeptide chains** stored in a directory:
```
python3 ./bin/repeatsdb-lite-2.py -e 1 -i /input/directory -o /output/directory -t /temporary/directory
```

If you want to download a specific structure from PDB (e.g. chain C of 4g8l PDB structure)
```
python3 ./bin/repeatsdb-lite-2.py -e 0 --pdb_id 4g8l --pdb_chain C -o /output/directory -t /temporary/directory
```

If you want to download a predicted structure from AlphaFold (e.g. UniProt ID: Q05823)
```
python3 ./bin/repeatsdb-lite-2.py -e 0 --uniprot_id Q05823 --af_version 4 -o /output/directory -t /temporary/directory
```
