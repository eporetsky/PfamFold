import pandas as pd

def load_pfam_data(filepath):
    """
    Loads the Pfam data from a TSV file.
    
    Args:
        filepath (str): Path to the Pfam TSV file.
    
    Returns:
        pd.DataFrame: DataFrame containing Pfam data.
    """
    return pd.read_csv(filepath, sep="\t")

def filter_pfam_domains(pfam_all, min_genes=20):
    """
    Filters Pfam domains based on the number of associated genes.
    
    Args:
        pfam_all (pd.DataFrame): DataFrame containing Pfam data.
        min_genes (int): Minimum number of genes for a domain to be included.
    
    Returns:
        pd.DataFrame: Filtered DataFrame with domain counts and names.
    """
    pfam_counts = pfam_all.groupby("domain_id").count()
    pfam_counts = pfam_counts.sort_values("gene", ascending=False).reset_index()
    pfam_counts = pfam_counts[["domain_id", "gene"]]
    pfam_names = pfam_all[["domain_id", "domain_name"]].drop_duplicates()
    pfam_counts = pfam_counts.merge(pfam_names, on="domain_id", how="left")
    return pfam_counts[pfam_counts["gene"] > min_genes]
