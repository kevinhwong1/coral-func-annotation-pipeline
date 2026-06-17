# Quick-start: load the coral annotation table in R.
library(tidyverse)

# Load annotation table (protein_sequence excluded from GitHub-hosted file)
ann <- read_tsv("species/Gfas/annotation/Gfas_master_annotation.tsv")
cat("Loaded:", nrow(ann), "genes x", ncol(ann), "columns\n")

# Annotation status breakdown
ann %>% count(annotation_status)

# Confidence tier breakdown
ann %>% count(confidence_tier)

# Gold-tier secreted ligands
ligands <- ann %>% filter(confidence_tier == "Gold", cellchat_role == "Ligand")
cat("Gold-tier ligands:", nrow(ligands), "\n")

# Parse pipe-delimited GO terms into long format
go_long <- ann %>%
  select(gfas_gene_id, all_go_terms) %>%
  separate_rows(all_go_terms, sep = "\\|") %>%
  filter(!is.na(all_go_terms))

# Join to Seurat metadata (example)
# seurat_meta <- FetchData(seurat_obj, vars = c("gene_id", "seurat_clusters"))
# merged <- seurat_meta %>% left_join(ann, by = c("gene_id" = "gfas_gene_id"))
