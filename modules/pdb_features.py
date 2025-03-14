# Add classify_structure and classify_all_structures functions for secondary structure classification
import os
import glob
from Bio import PDB

def calculate_average_plddt(pdb_file):
    """
    Calculate the average pLDDT score from a PDB file.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        float: The average pLDDT score.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    plddt_scores = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    plddt_scores.append(atom.bfactor)

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    return None

def classify_stride(pdb_file, threshold=0.1):
    """
    Classify a protein structure as 'none', 'a-helix only', 'b-sheet only', or 'mixed'
    based on secondary structure content using STRIDE results.
    
    Parameters:
    - pdb_file: str, the filename of the PDB file (e.g., "species_proteinid_pfamid_start-end.pdb").
    - threshold: float, the minimum proportion of residues required to classify as a specific structure.
    
    Returns:
    - str, the classification of the domain.
    """
    # Parse the filename to extract species, protein ID, Pfam domain, and residue range
    filename = os.path.basename(pdb_file)
    parts = filename.split("_")
    if len(parts) < 3:
        raise ValueError(f"Invalid filename format: {filename}")
    
    species = parts[0]
    protein_id = parts[1]
    residue_range = parts[3].replace(".pdb", "")
    start_res, end_res = map(int, residue_range.split("-"))

    # Locate the corresponding STRIDE result file
    stride_file = f"input/stride/{species}/{protein_id}.stride"

    # Parse the STRIDE result file
    helix_count = 0
    sheet_count = 0
    total_count = 0
    with open(stride_file, "r") as f:
        for line in f:
            if line.startswith("ASG"):  # Look for residue-specific secondary structure assignments
                parts = line.split()
                res_num = int(parts[3])  # Residue number
                structure = parts[5]    # Secondary structure type (e.g., 'H', 'E', 'C', etc.)

                # Check if the residue is within the specified range
                if start_res <= res_num <= end_res:
                    total_count += 1
                    if structure in ('H', 'G', 'I'):  # Helices
                        helix_count += 1
                    elif structure in ('E', 'B'):    # Beta-sheets
                        sheet_count += 1

    # Avoid division by zero if no residues are found in the range
    if total_count == 0:
        return "none"

    # Calculate proportions
    helix_proportion = helix_count / total_count
    sheet_proportion = sheet_count / total_count
    other_proportion = (total_count - helix_count - sheet_count) / total_count

    # Classify the structure based on the threshold
    if helix_proportion >= threshold and sheet_proportion < threshold:
        seconadry_category = "a-helix"
    elif sheet_proportion >= threshold and helix_proportion < threshold:
        seconadry_category = "b-sheet"
    elif helix_proportion >= threshold and sheet_proportion >= threshold:
        seconadry_category = "mixed"
    else:
        seconadry_category = "coil"
    
    return(seconadry_category, helix_proportion, sheet_proportion, other_proportion)

def process_pdb_features(input_path, output_path, num_threads):
    """
    Process PDB files from a single file or a folder in parallel and save the average pLDDT scores to a TSV file.

    Parameters:
        input_path (str): Path to tmp selres directory with PDBs
        output_path (str): Path to the output TSV file.
        num_threads (int): Number of threads to use for parallel processing.
    """
    pdb_files = [fl for fl in glob.glob(f"{input_path}/*.pdb")]

    pdb_features = []

    for pdb_file in pdb_files:
        avg_plddt = calculate_average_plddt(pdb_file)
        stride, helix, sheet, coil = classify_stride(pdb_file)
        pdb_features.append((pdb_file, avg_plddt, stride, helix, sheet, coil))

    with open(output_path, 'w') as out_file:
        out_file.write("Filename\tAverage_pLDDT\tSTRIDE\thelix\tsheet\tcoil\n")
        for pdb_file, avg_plddt, stride, helix, sheet, coil in pdb_features:
            out_file.write(f"{pdb_file}\t{avg_plddt:.2f}\t{stride}\t{helix}\t{sheet}\t{coil}\n")