import os
import shutil

def create_directories(output_path):
    """
    Creates the necessary directory structure for the project.
    
    Args:
        project_name (str): Name of the project directory.
    """
    os.makedirs("tmp", exist_ok=True)
    os.makedirs(output_path, exist_ok=True)
    subdirs = ["cluster_reps_foldseek", 
               "cluster_scores_foldseek",
               "cluster_reps_agg",
               "cluster_scores_agg",
               "plots", 
               "selres", 
               "pdb_features", 
               "qc_groups"]
    for subdir in subdirs:
        os.makedirs(f"{output_path}/{subdir}", exist_ok=True)

def clean_tmp_directory(pfam_id):
    """
    Cleans the temporary directory.
    """
    shutil.rmtree(f"tmp/{pfam_id}", ignore_errors=True)
    os.makedirs(f"tmp/{pfam_id}", exist_ok=True)