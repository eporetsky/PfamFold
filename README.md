# PfamFold

PfamFold is a Python-based workflow for analyzing structural variability within protein domains. By leveraging predicted protein structures from AlphaFold2 and domain annotations from the Pfam database, this tool identifies and clusters structural variations within Pfam families. The workflow incorporates tools like FoldSeek and agglomerative clustering to explore structural variability, providing insights into domain curation and refinement.

## Features

- Extracts and processes protein domain structures based on Pfam annotations.
- Computes structural features such as pLDDT scores and secondary structure classifications.
- Performs clustering of domain structures using FoldSeek and agglomerative clustering.
- Generates TM-score matrices and outputs clustered domain representatives for further analysis.

## Requirements

The dependencies for this project are managed via Conda. The required environment can be set up using the included `environment.yml` file.

### Setting up the environment

Run the following commands to create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate pfamfold
```

## Usage

The main script is `main.py`, which processes Pfam domains and outputs results in the `results/` directory. 

### Steps:

1. Place the input files (e.g., Pfam annotations, PDB files) in the appropriate directories as required by the workflow.
2. Run the script using:

```bash
python main.py
```

3. Results, including structural features, clustering outputs, and representative domains, will be stored in the `results/` directory.

### Input Requirements:

The workflow requires an `input/` directory structured as follows:
- **input/pdb/**: Contains the predicted AlphaFold2 structures.
- **input/selres/**: Contains the structural portion of the Pfam domains.
- **input/stride/**: Contains the secondary structure annotation files.
- **input/interpro.pfam.tsv***: Pfam domain data in TSV format.

Each of the directories must include subdirectories for each genome being analyzed. For example, if analyzing multiple genomes, each genome should have its own subdirectory within `pdb/`, `selres/`, and `stride/`, containing the relevant files for that genome.

### Output:
- Structural features (e.g., pLDDT scores, secondary structure classifications).
- Clustered domain representatives.
- TM-score matrices for structural comparisons.

## Key Tools and Libraries

- **AlphaFold2**: Predicted protein structures.
- **Pfam**: Domain annotations.
- **FoldSeek**: Clustering based on structural similarity.
- **pandas** & **scikit-learn**: Data analysis and clustering.
- **Biopython**: Handling PDB files.
- **STRIDE**: Calculating secondary structure

## Contributing

Contributions and suggestions are welcome. Feel free to submit issues or pull requests to improve the project.

## License

This project is licensed under the MIT License. See the LICENSE file for details.