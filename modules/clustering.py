import os
import shutil
import pandas as pd
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage

def run_foldseek(input_dir, output_dir, tmp_dir, pfam_id, clustering_threshold=0.5, threads=16):
    """
    Runs Foldseek clustering on Pfam PDB files.
    
    Args:
        input_dir (str): Directory containing input PDB files.
        output_dir (str): Directory to save Foldseek output.
        tmp_dir (str): Temporary directory for Foldseek.
        pfam_id (str): Pfam domain ID.
        clustering_threshold (float): Clustering threshold for Foldseek.
        threads (int): Number of threads to use.
    """
    os.system(f'foldseek easy-cluster {input_dir} {output_dir}/{pfam_id} {tmp_dir} -c {clustering_threshold} --threads {threads} > /dev/null')

def generate_distance_matrix(foldseek_output, output_file):
    """
    Generates a distance matrix from Foldseek output.
    
    Args:
        foldseek_output (str): Path to the Foldseek output file.
        output_file (str): Path to save the distance matrix.
    """
    foldseek = pd.read_csv(foldseek_output, sep="\t", header=None)
    foldseek = foldseek[[0, 1, 3]]  # Extract relevant columns
    foldseek['distance'] = foldseek[3]
    wide_df = foldseek.pivot(index=0, columns=1, values='distance').fillna(1)
    wide_df.to_csv(output_file, sep="\t", index=False)

def generate_tm_score_matrix(pfam_id, rep_dir, all_dir, output_file, threads=16):
    """
    Generate a TM-score matrix using Foldseek easy-search between cluster representatives and all members.

    Args:
        rep_dir (str): Directory containing cluster representative PDB files.
        all_dir (str): Directory containing all PDB files.
        output_file (str): Path to save the TM-score matrix in wide format.
        threads (int): Number of threads to use for Foldseek.
    """

    # Use the representatives of the Foldseek clusters to calculate RMSDs
    # based on reps x all. This will save on resources for large sets
    clusters_file = output_file.replace("cluster_scores_foldseek","cluster_reps_foldseek").replace(".wide.tsv", "_cluster.tsv")
    foldseek = pd.read_csv(f"{clusters_file}", sep="\t", header=None)
    grouped = foldseek.groupby(0).count()
    singletons = grouped[grouped[1]==1]
    if len(singletons) > 1000:
        grouped = grouped[grouped[1]>1]
    rep_groups = list(grouped.index)
    if len(rep_groups) < 2:
        return(True)
    for rep in rep_groups:
        shutil.copy(f"tmp/{pfam_id}/selres/{rep}.pdb", f"tmp/{pfam_id}/selres_clust/")

    # Temporary output file for Foldseek results
    tmp_output = output_file.replace(".wide.tsv", ".aln.tsv")

    # Run Foldseek easy-search
    os.system(
        f'foldseek easy-search {rep_dir} {all_dir} {tmp_output} tmp/{pfam_id} '
        f'--alignment-type 1 --tmscore-threshold 0.0 --format-output query,target,rmsd,alntmscore,u,t '
        f'--exhaustive-search 1 -e inf --threads {threads} > /dev/null'
    )

    # Read the Foldseek output
    foldseek = pd.read_csv(tmp_output, sep="\t", header=None)
    foldseek = foldseek[[0, 1, 3]]  # Extract query, target, and TM-score columns
    foldseek.columns = ["query", "target", "tm_score"]

    # Convert similarity to distance (optional, here we keep TM-scores as is)
    foldseek["distance"] = foldseek["tm_score"]

    # Create a wide-format distance matrix
    wide_df = foldseek.pivot(index="query", columns="target", values="distance")
    wide_df = wide_df.fillna(1)  # Fill NaN with 1 (distance of 1 for missing pairs)
    wide_df.to_csv(output_file, sep="\t")

    # Clean up temporary files
    #os.remove(tmp_output)

def perform_agglomerative_clustering(tm_score_file, 
                                     pdb_features_file,
                                     input_foldseek_cluster,
                                     output_path_agg_reps,
                                     output_path_agg_tm, 
                                     distance_threshold=0.4):
    """
    Perform agglomerative clustering on TM-score data and select representatives based on pLDDT scores.

    Inputs:
        tm_score_file (str): Path to the TM-score matrix file.
        pdb_features_file (str): Path to the pLDDT scores file.
        input_foldseek_cluster (str): Path to the FoldSeek cluster file.
        output_path_agg_reps (str): Path to save the agglomerative clustering reps.
        output_path (str): Path to save the agglomerative clustering TM-scores.
        distance_threshold (float): Distance threshold for clustering (default: 0.4, which is 1 - TM-score).
    """
    # Read the TM-score data
    tm_df = pd.read_csv(tm_score_file, sep="\t", index_col=0).T
    tm_df = tm_df.loc[tm_df.columns, tm_df.columns]  # Ensure symmetry
    # Read the pLDDT scores, omitting other columns
    plddt = pd.read_csv(pdb_features_file, sep="\t")[['Filename','Average_pLDDT']]
    plddt_data = plddt.set_index('Filename')['Average_pLDDT'].to_dict()

    # Convert TM-scores to distances (1 - TM-score)
    distance_matrix = 1 - tm_df.values

    # Perform Agglomerative Clustering
    clustering = AgglomerativeClustering(
        linkage='single',
        metric='precomputed',
        distance_threshold=distance_threshold,
        n_clusters=None
    )
    cluster_labels = clustering.fit_predict(distance_matrix)

    # Map cluster labels to representatives
    cluster_to_representatives = {}
    for idx, cluster_id in enumerate(cluster_labels):
        representative = tm_df.columns[idx]
        if cluster_id not in cluster_to_representatives:
            cluster_to_representatives[cluster_id] = []
        cluster_to_representatives[cluster_id].append(representative)

    # Select the representative with the highest pLDDT for each cluster
    selected_representatives = {}
    reduced_mapping_dict = {}
    for cluster_id, representatives in cluster_to_representatives.items():
        best_rep = max(representatives, key=lambda rep: plddt_data.get(rep, 0))
        selected_representatives[cluster_id] = best_rep
        for clust_member in representatives:
            reduced_mapping_dict[clust_member] = best_rep

    # Generate the agglomerative clustering cluster representatives and members table
    # Replace all values in the DataFrame using the dictionary
    reduced_clusters = pd.read_csv(input_foldseek_cluster, sep="\t", header=None)
    reduced_clusters[0] = reduced_clusters[0].replace(reduced_mapping_dict)
    reduced_clusters.columns = ["representatives", "members"]
    reduced_clusters.to_csv(output_path_agg_reps, sep="\t", index=None, header=None)

    # Save the reduced clusters
    reduced_tm_df = tm_df.loc[selected_representatives.values(), selected_representatives.values()]
    reduced_tm_df.to_csv(output_path_agg_tm, sep="\t")

    return reduced_tm_df