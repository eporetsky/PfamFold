import pandas as pd
import os

fl_df = pd.read_csv("interpro.dedup.pfam.tsv", sep="\t")
fl_df = fl_df.drop_duplicates("gene", axis=1)
fl_df = fl_df[["genome", "gene"]]

output_df = []
for index, row in fl_df.iterrows():
    genome = row["genome"]
    gene = row['gene']
    input_path = os.path.join("pdb", genome, f"AF-{gene}-F1-model_v4.pdb")
    output_path = os.path.join("stride", genome, f"{gene}.stride")
    output_df.append([input_path, output_path])
output_df = pd.DataFrame(output_df)
output_df.to_csv("stride.tsv", sep="\t", index=None, header=None)

os.makedirs(os.path.join("stride"), exist_ok=True)
for genome in fl_df["genome"].unique().tolist():
    os.makedirs(os.path.join("stride", genome), exist_ok=True)
    
os.system("parallel -j70 --colsep '\t' -a stride.tsv './stride {1} > {2}'")
