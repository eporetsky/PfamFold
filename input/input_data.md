# Input Data Preparation Guide

This document provides instructions on how to download and prepare the input data required for the workflow. The data includes AlphaFold2-predicted protein structures, which are processed and analyzed in the pipeline.

## Conda environment

The dependencies for this project are managed via Conda. Run the following commands to create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate pfamfold
```

### Getting the AlphaFold2 PDB files

Below are the scripts to recreate the input data used in the analysis. The AlphaFold2 PDB structures for the 16 reference genomes can be downloaded from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/download). 

Ensure that the `input/` directory is structured as follows:

```
input/
├── pdb/
│   ├── celegans/
│   ├── calbicans/
│   ├── zebrafish/
│   └── ... (other genomes)
├── fasta/
│   ├── celegans.fasta
│   ├── calbicans.fasta
│   ├── zebrafish.fasta
│   └── ... (other genomes)
└── interproscan/
    ├── celegans.tsv
    ├── calbicans.tsv
    ├── zebrafish.tsv
    └── ... (other genomes)
```

## Steps to Prepare Input Data

### 1. Organize the Files

Create a directory structure to store the downloaded files:

```bash
mkdir pdb
cd pdb

# Create subdirectories for each genome
mkdir celegans calbicans zebrafish ...
```

### 2. Download AlphaFold2 PDB Files

Use the following commands to download the AlphaFold2-predicted protein structures for various genomes:

```bash
axel -n 32 https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000001940_6239_CAEEL_v4.tar
axel -n 32 https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000559_237561_CANAL_v4.tar
axel -n 32 https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000437_7955_DANRE_v4.tar
# (Repeat for all other genomes listed in the script)
```

### 3. Extract the Files

Extract the downloaded .tar files into their respective directories:

```bash
tar -xf UP000001940_6239_CAEEL_v4.tar -C celegans
tar -xf UP000000559_237561_CANAL_v4.tar -C calbicans
tar -xf UP000000437_7955_DANRE_v4.tar -C zebrafish
# (Repeat for all other genomes)
```

### 4. Remove Unnecessary Files

Remove .cif.gz files, keeping only the .pdb.gz files:

```bash
rm celegans/*.cif.gz
rm calbicans/*.cif.gz
rm zebrafish/*.cif.gz
# (Repeat for all other directories)
```

### 5. Decompress the PDB Files

Decompress the .pdb.gz files using pigz:

```bash
pigz -d celegans/*
pigz -d calbicans/*
pigz -d zebrafish/*
# (Repeat for all other directories)
```

### 6. Generate FASTA Files

Run the pdb2fasta.db.py script to generate a single FASTA file for each genome:

```bash
# Go back to the input/ directory
cd ..
mkdir fasta

python pdb2fasta.db.py pdb/celegans/ fasta/celegans.fasta
python pdb2fasta.db.py pdb/calbicans/ fasta/calbicans.fasta
python pdb2fasta.db.py pdb/zebrafish/ fasta/zebrafish.fasta
# (Repeat for all other genomes)
```

### 7. Run InterProScan

Use InterProScan to analyze the FASTA files. Replace interproscan-5.68-100.0-64 with the path to your InterProScan installation:

```bash
mkdir interproscan

bash interproscan-5.68-100.0-64/interproscan.sh -i fasta/celegans.fasta --output-dir interproscan --formats TSV --disable-residue-annot --applications Pfam
bash interproscan-5.68-100.0-64/interproscan.sh -i fasta/calbicans.fasta --output-dir interproscan --formats TSV --disable-residue-annot --applications Pfam
bash interproscan-5.68-100.0-64/interproscan.sh -i fasta/zebrafish.fasta --output-dir interproscan --formats TSV --disable-residue-annot --applications Pfam
# (Repeat for all other genomes)
```

Note: I'd recommend using the latest version of [InterProScan](https://www.ebi.ac.uk/interpro/download/interproscan/) to be able to use the online lookup service to speed up the analysis. Running an older version (such as 5.68-100.0-64) runs all predictions on your local machine. Also recommended to adjust the memory and threads when running InterProScan.

### 8. Combine InterProScan results

Use the `input_ips.py` python script to combine all the InterProScan results for the separate genomes into a single file:

```bash
python input_ips.py
```

Note: This script assumes gene names are all given in the AlphaFold2 nomenclature (AF-P18843-F1-model_v4) and will extract the UniProt portion of it (P18843). It will also filter out all UniProt IDs with more than one predicted AlphaFold2 structure.

### 9. Extract the Pfam domain structures

Use the `input_selres.py` python script to read the `interpro.dedup.pfam.tsv` file and the selres function from PDBtools to extract the Pfam domain structures from the complete protein structure.

```bash
python input_selres.py
```


### 10. Get the STRIDE secondary structures

Use the `input_stride.py` python script to assign secondary structures to the full AlphaFold2 PDB files.

```bash
python input_stride.py
```

Note: The executable file for STRIDE is not included but can be downloaded from [their website](https://ssbio.readthedocs.io/en/latest/instructions/stride.html).