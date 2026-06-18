# Column Descriptions

## Contents

- [Gene Identity](#gene-identity)
- [eggNOG-mapper](#eggnog-mapper)
- [OrthoFinder](#orthofinder)
- [InterProScan](#interproscan)
- [Reciprocal Best Hit (RBH) vs SwissProt](#reciprocal-best-hit-rbh-vs-swissprot)
- [SignalP](#signalp)
- [Tier 3: Custom HMM Profiles](#tier-3-custom-hmm-profiles)
- [Tier 4: Curated Coral Sequences](#tier-4-curated-coral-sequences)
- [Genomic Coordinates](#genomic-coordinates)
- [DeepTMHMM](#deeptmhmm)
- [Protein Properties](#protein-properties)
- [Classification](#classification)
- [Derived Annotations](#derived-annotations)
- [Delimiter Summary](#delimiter-summary)

---

Master annotation table schema. Standard species (Gfas, Pdam, Amil, Spis, Opat, Xspp, Nvec): 53 columns. NCBI-sourced species (Acer, Ofav, Amur): 54 columns, with `{species}_locus_tag` added at column 2.

---

## Gene Identity

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 1 | `{species}_gene_id` | string | Primary key; gene model ID, species-prefixed (e.g. `gfas_gene_id`, `pdam_gene_id`) | `gfas1.m1.1.m1` |

---

## eggNOG-mapper

Functional annotations from **eggNOG-mapper v2** ([GitHub](https://github.com/eggnogdb/eggnog-mapper)) against the eggNOG 5.0 database.

> Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. *Mol Biol Evol.* 2021;38(12):5825–5829. https://doi.org/10.1093/molbev/msab293

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 2 | `eggnog_ortholog` | string | Best eggNOG seed ortholog match (taxon.proteinID) | `45351.EDO49519` |
| 3 | `eggnog_human_gene` | string | Preferred human gene name from the matched OG | `WDR76` |
| 4 | `eggnog_description` | string | Functional description from eggNOG OG | `regulation of DNA damage checkpoint` |
| 5 | `eggnog_go_terms` | string | GO terms from eggNOG (comma-delimited) | `GO:0003674,GO:0005634` |
| 6 | `eggnog_kegg_ko` | string | KEGG KO number(s) (comma-delimited) | `ko:K15362` |
| 7 | `eggnog_kegg_pathway` | string | KEGG pathway IDs (comma-delimited) | `ko03440,ko03460` |
| 8 | `eggnog_cog_cat` | string | COG functional category letter code | `S` |
| 9 | `eggnog_pfams` | string | Pfam domain names from eggNOG (comma-delimited) | `WD40` |

---

## OrthoFinder

Ortholog assignments from **OrthoFinder v2** ([GitHub](https://github.com/davidemms/OrthoFinder)) run in MSA mode across 4 species (coral, human, mouse, *Nematostella vectensis*).

> Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biol.* 2019;20:238. https://doi.org/10.1186/s13059-019-1832-y

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 10 | `of_orthogroup` | string | OrthoFinder orthogroup ID | `OG0006633` |
| 11 | `of_human_ortholog` | string | Human ortholog UniProt entry name | `WDR76_HUMAN` |
| 12 | `of_mouse_ortholog` | string | Mouse ortholog UniProt entry name | `WDR76_MOUSE` |
| 13 | `of_nvec_ortholog` | string | *Nematostella vectensis* ortholog ID | `NV2t012202001.1` |
| 14 | `of_relationship` | string | Orthology relationship type | `one-to-one` / `one-to-many` / `many-to-many` |

---

## InterProScan

Domain and family annotations from **InterProScan 5.78-109.0** ([EBI](https://www.ebi.ac.uk/interpro/)) across 17 member databases.

> Jones P, Binns D, Chang HY, et al. InterProScan 5: genome-scale protein function classification. *Bioinformatics.* 2014;30(9):1236–1240. https://doi.org/10.1093/bioinformatics/btu031

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 15 | `interpro_pfam` | string | Pfam accession(s) from IPS (semicolon-delimited) | `PF00400` |
| 16 | `interpro_panther` | string | PANTHER accession(s) from IPS (semicolon-delimited) | `PTHR14773` |
| 17 | `interpro_go_terms` | string | GO terms from IPS (comma-delimited) | `GO:0003677,GO:0005515` |

---

## Reciprocal Best Hit (RBH) vs SwissProt

Reciprocal best BLAST hits against UniProt Swiss-Prot using **DIAMOND v2.2.1** ([GitHub](https://github.com/bbuchfink/diamond)) (blastp, e-value ≤ 1e-5).

> Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods.* 2021;18:366–368. https://doi.org/10.1038/s41592-021-01101-x

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 18 | `rbh_swissprot_acc` | string | SwissProt accession of RBH | `Q4KLQ5` |
| 19 | `rbh_gene_name` | string | Gene name of RBH (GENE_SPECIES format) | `WDR76_XENLA` |
| 20 | `rbh_pct_identity` | float | Percent identity to RBH | `38.89` |
| 21 | `rbh_evalue` | float | E-value of RBH alignment | `3.13e-49` |

---

## SignalP

Signal peptide predictions from **SignalP 6.0** ([DTU server](https://services.healthtech.dtu.dk/services/SignalP-6.0/)) in slow-sequential mode.

> Teufel F, Almagro Armenteros JJ, Johansen AR, et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat Biotechnol.* 2022;40:1023–1025. https://doi.org/10.1038/s41587-021-01156-3

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 22 | `signalp_prediction` | string | Signal peptide class | `SP` / `OTHER` / `LIPO` / `TAT` |
| 23 | `signalp_prob` | float | Signal peptide probability (0–1) | `0.97` |

---

## Tier 3: Custom HMM Profiles

Hits against a curated set of 179 coral/metazoan functional HMM profiles using **HMMER 3.4 hmmsearch** ([hmmer.org](http://hmmer.org/)) covering AMPs, venom/toxin, immunity, transcription factors, and ECM domains.

> Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol.* 2011;7(10):e1002195. https://doi.org/10.1371/journal.pcbi.1002195

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 24 | `tier3_domain_hits` | string | Matched HMM profile name(s) (comma-delimited) | `Pentaxin` |
| 25 | `tier3_best_evalue` | float | Best e-value across all Tier 3 profile hits | `6.1e-33` |

---

## Tier 4: Curated Coral Sequences

BLAST hits against a curated set of 17 coral-specific reference sequences (SCRiPs, SOMPs, neuropeptides, cnidocyte proteins) using **DIAMOND v2.2.1**.

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 26 | `tier4_blast_hit` | string | Reference sequence ID of best Tier 4 hit | `SCRiP_Venom_Amil_SCRiP1_B7X8E0` |
| 27 | `tier4_source` | string | Tier 4 family | `SCRiP` / `SOMP` / `Neuropep` / `Cnidocyte` |
| 28 | `tier4_pct_identity` | float | Percent identity to Tier 4 hit | `31.4` |
| 29 | `tier4_evalue` | float | E-value of Tier 4 BLAST hit | `4.1e-08` |
| 30 | `tier4_bitscore` | float | Bitscore of Tier 4 BLAST hit | `44.3` |

---

## Genomic Coordinates

Parsed from the species reference GFF3 file.

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 31 | `scaffold` | string | Scaffold or chromosome ID | `Sc0000000` |
| 32 | `gene_start` | int | Gene start position (1-based) | `28422` |
| 33 | `gene_end` | int | Gene end position (1-based) | `41700` |
| 34 | `strand` | string | Strand orientation | `+` / `-` |
| 35 | `gene_length_bp` | int | Gene span in base pairs (end − start) | `13279` |
| 36 | `n_exons` | int | Number of exons from GFF3 | `13` |

---

## DeepTMHMM

Transmembrane topology predictions from **DeepTMHMM v1.0** ([DTU server](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/)). Run as standalone academic license. Proteomes split into 2,000-sequence chunks.

> Hallgren J, Tsirigos KD, Pedersen MD, et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. *bioRxiv.* 2022. https://doi.org/10.1101/2022.04.08.487609

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 37 | `tmhmm_topology` | string | Predicted topology class | `cytoplasmic` / `TMhelix` / `signal+TMhelix` |
| 38 | `tmhmm_n_tmrs` | int | Number of predicted TM helices | `0` / `7` |

---

## Protein Properties

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 39 | `protein_length_aa` | int | Protein length in amino acids | `512` |
| 40 | `protein_sequence` | string | Full amino acid sequence | `MPRSDTG...` |

---

## Classification

Functional classification assigned by the pipeline merge script (`08_merge_annotate.py`). See [classification_system.md](classification_system.md) for full logic.

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 41 | `confidence_tier` | string | Evidence-based confidence tier | `Gold` / `Silver` / `Bronze` / `Informative` / `Unclassified` |
| 42 | `protein_group` | string | Primary functional group (24 groups; for Seurat/plotting) | `GPCR` |
| 43 | `protein_subgroup` | string | Specific subclass within protein_group | `Pentraxin_CRP` |
| 44 | `protein_group_secondary` | string | Secondary functional group for dual-role genes (NaN if single-role) | `Neuron_Neuropeptide` |
| 45 | `cellchat_role` | string | Predicted intercellular communication role | `Ligand` / `Receptor` / `Receptor_candidate` / `Neither` |

---

## Derived Annotations

Computed or merged fields added during the final annotation step using **BioPython v1.87** ([biopython.org](https://biopython.org/)) and the **KEGG REST API** (https://rest.kegg.jp/).

> Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 46 | `best_human_gene_name` | string | Best human gene symbol (OrthoFinder > RBH > eggNOG, in priority order) | `WDR76` |
| 47 | `all_go_terms` | string | Merged, deduplicated GO terms from all sources (**pipe-delimited**) | `GO:0003674\|GO:0005634` |
| 48 | `cog_description` | string | Full text description of COG category | `Function unknown` |
| 49 | `interpro_pfam_names` | string | Human-readable Pfam domain names from IPS (semicolon-delimited) | `WD domain, G-beta repeat` |
| 50 | `annotation_status` | string | Annotation completeness per Stephens et al. 2026 GBE | `Light` / `Partial` / `Dark` |
| 51 | `isoelectric_point` | float | Theoretical isoelectric point computed by BioPython ProteinAnalysis | `8.48` |
| 52 | `molecular_weight_kda` | float | Theoretical molecular weight in kDa computed by BioPython ProteinAnalysis | `57.31` |
| 53 | `kegg_pathway_names` | string | Human-readable KEGG pathway names fetched via KEGG REST API (semicolon-delimited) | `Homologous recombination` |

---

## Delimiter Summary

| Column(s) | Delimiter |
|-----------|-----------|
| `eggnog_go_terms`, `interpro_go_terms` | comma `,` |
| `eggnog_kegg_ko`, `eggnog_kegg_pathway`, `eggnog_pfams`, `tier3_domain_hits` | comma `,` |
| `interpro_pfam`, `interpro_panther`, `interpro_pfam_names`, `kegg_pathway_names` | semicolon `;` |
| `all_go_terms` | pipe `\|` |

> **Note:** `all_go_terms` uses pipe-delimited values because the source GO columns use commas; pipes avoid ambiguity when parsing the merged field.
