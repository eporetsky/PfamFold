import os
import pandas as pd
import shutil
from modules.file_utils import create_directories, clean_tmp_directory
from modules.pfam_analysis import load_pfam_data, filter_pfam_domains
from modules.clustering import run_foldseek, generate_tm_score_matrix, perform_agglomerative_clustering
from modules.pdb_features import process_pdb_features

# Process each Pfam domain
def pfamfold(pfam_id):
    global pfam_all
    #print(f"Processing Pfam domain: {pfam_id}")
    output_path = "results"
    try:
        # Create temporary directories for processing
        clean_tmp_directory("pfam_id")
        os.makedirs(f"tmp/{pfam_id}/selres", exist_ok=True)
        os.makedirs(f"tmp/{pfam_id}/selres_clust", exist_ok=True)

        # Extract Pfam domain-specific PDB files
        tmp_df = pfam_all[pfam_all["domain_id"] == pfam_id].copy()
        count = 0

        for _, gene_row in tmp_df.iterrows():
            genotype = gene_row["genotype"]
            gene = gene_row["gene"]
            start = gene_row["start"]
            end = gene_row["end"]

            input_path = f"input/selres/{genotype}/{genotype}_{gene}_{pfam_id}_{start}-{end}.pdb"
            output_file = f"tmp/{pfam_id}/selres"
            os.system(f"cp {input_path} {output_file}")
            count += 1

        # Calculate pLDDT scores
        process_pdb_features(
            input_path=f"tmp/{pfam_id}/selres",
            output_path=f"{output_path}/pdb_features/{pfam_id}.tsv",
            num_threads=os.cpu_count()
        )

        # Run Foldseek clustering
        run_foldseek(
            input_dir=f"tmp/{pfam_id}/selres/",
            output_dir=f"{output_path}/cluster_reps_foldseek",
            tmp_dir=f"tmp/{pfam_id}",
            pfam_id=pfam_id
        )

        # Generate TM-score matrix
        tm_score_matrix_file = f"{output_path}/cluster_scores_foldseek/{pfam_id}.wide.tsv"
        no_clusters = generate_tm_score_matrix(
            pfam_id=pfam_id,
            rep_dir=f"tmp/{pfam_id}/selres_clust",
            all_dir=f"tmp/{pfam_id}/selres",
            output_file=tm_score_matrix_file,
            threads=8 #os.cpu_count()
        )

        shutil.rmtree(f"tmp/{pfam_id}", ignore_errors=True)
        #clean_tmp_directory("pfam_id")

        if no_clusters:
            return(None)

        # Perform agglomerative clustering
        reduced_tm_df = perform_agglomerative_clustering(
            tm_score_file=tm_score_matrix_file,
            pdb_features_file=f"{output_path}/pdb_features/{pfam_id}.tsv",
            input_foldseek_cluster=f"{output_path}/cluster_reps_foldseek/{pfam_id}_cluster.tsv",
            output_path_agg_reps=f"{output_path}/cluster_reps_agg/{pfam_id}_cluster.tsv",
            output_path_agg_tm=f"{output_path}/cluster_scores_agg/{pfam_id}.tsv"
        )
        if len(reduced_tm_df) == 1:
            return(None)

    except Exception as e:
        # Handle the exception
        print(f"{pfam_id}: {e}")

if __name__ == "__main__":
    os.system("rm -r results/")
    os.system("rm -r tmp/*")

    # Initialize output path
    output_path = "results"
    create_directories(output_path)

    # Load and filter Pfam data
    pfam_all = load_pfam_data("input/interpro.pfam.tsv")
        
    # Filters Pfam domains
    pfam_counts = filter_pfam_domains(pfam_all, min_genes=20)
    
    data = pfam_counts["domain_id"].tolist()

    from multiprocessing import Pool
    # Set the maximum number of processes
    max_processes = 24

    with Pool(max_processes) as pool:
        pool.map(pfamfold, data)