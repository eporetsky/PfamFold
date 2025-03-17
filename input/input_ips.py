import os
import glob
import pandas as pd

ips_all = []
for fl in glob.glob("interproscan/*.fa.tsv"):
    # Load the data
    ips = pd.read_csv(fl, sep='\t', header=None)
    genotype = os.path.basename(fl).replace(".fa.tsv", "")
    # Rename columns based on the provided example
    columns = ['af2_id', 'transcript_id', 'length', 'database', 'domain_id', 'domain_name', 
            'start', 'end', 'evalue', 'status', 'date', 'entry_type', 'entry_name', "skip1", "skip2"]
    ips.columns = columns
    # Filter Pfam annotations
    ips = ips[ips['database'] == 'Pfam']

    # Add the genotype and gene columns
    ips['genotype'] = genotype
    ips["gene"] = ips["af2_id"].str.split("-").str[1] # AF-P18843-F1-model_v4

    # Find and remove all gene IDs that have more than one unique AF2 structure
    ips_unique = ips["af2_id"].unique().tolist()
    ips_unique = [gene_id.split("-")[1] for gene_id in ips_unique]
    duplicates = [gene for gene in ips_unique if ips_unique.count(gene) > 1]
    ips = ips[~ips["gene"].isin(duplicates)]

    # Reorder columns and save table
    ips = ips[["genotype","gene","length","domain_id","domain_name","start","end","evalue"]]
    ips_all.append(ips)

# Concatenate all tables and save into a single file
ips_all = pd.concat(ips_all)
ips_all.to_csv("interpro.dedup.pfam.tsv", sep="\t", index=False)