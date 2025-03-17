import os
import pandas as pd

pfam_all = pd.read_csv("interpro.dedup.pfam.tsv", sep="\t")

os.makedirs("selres", exist_ok=True)

output_df = []
for index, row in pfam_all.iterrows():
    genotype = row["genotype"]
    gene = row['gene']
    domain_id = row["domain_id"]
    start = row['start']
    end = row['end']
    input_path = os.path.join("pdb", genotype, f"AF-{gene}-F1-model_v4.pdb")
    output_path = os.path.join("selres", genotype, f"{genotype}_{gene}_{domain_id}_{start}-{end}.pdb")
    output_df.append([start, end, input_path, output_path])
output_df = pd.DataFrame(output_df)
output_df.to_csv("selres.pfam.tsv", sep="\t", index=None, header=None)

os.makedirs(os.path.join("large", "selres"), exist_ok=True)
for genotype in pfam_all["genotype"].unique().tolist():
    os.makedirs(os.path.join("large", "selres",genotype), exist_ok=True)
    
os.system("parallel -j32 --colsep '\t' -a selres.pfam.tsv 'pdb_selres -{1}:{2} {3} > {4}'")

os.system("rm selres.pfam.tsv")