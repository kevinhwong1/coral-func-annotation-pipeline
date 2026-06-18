"""
Quick-start: load the coral annotation table in Python.
"""
import pandas as pd

# Load annotation table (protein_sequence excluded from GitHub-hosted file)
df = pd.read_csv("species/Gfas/annotation/Gfas_master_annotation.tsv", sep="\t")
print(f"Loaded: {df.shape[0]} genes x {df.shape[1]} columns")

# Annotation status breakdown
print("\nAnnotation status:")
print(df["annotation_status"].value_counts())

# Confidence tier breakdown
print("\nConfidence tiers:")
print(df["confidence_tier"].value_counts())

# Gold-tier secreted ligands
ligands = df[(df["confidence_tier"] == "Gold") & (df["cellchat_role"] == "Ligand")]
print(f"\nGold-tier ligands: {len(ligands)}")

# Parse pipe-delimited GO terms into a list
df["go_list"] = df["all_go_terms"].str.split("|")

# Join to Seurat metadata (example)
# seurat_meta = pd.read_csv("seurat_metadata.csv", index_col=0)
# merged = seurat_meta.join(df.set_index("{species}_gene_id"), on="gene_id")  # replace {species}
