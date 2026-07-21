# *Xspp* Annotation Summary

## Contents

- [Run Details](#run-details)
- [Gene Universe](#gene-universe)
- [Confidence Tiers](#confidence-tiers)
- [Annotation Status](#annotation-status)
- [Protein Groups](#protein-groups)
- [CellChat Roles](#cellchat-roles)
- [Tool Versions Used](#tool-versions-used)
- [Gene Lists](#gene-lists)

---

## Run Details

| Field | Value |
|-------|-------|
| Genome version | GCF_021976095.1 |
| Proteome source | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_021976095.1/ |
| Proteome file | `xenia_protein.faa` |
| GFF3 file | `genomic.gff` |
| Pipeline run date | 2026-07-20 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Xspp_master_annotation.tsv` (26,116 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 26,116 |
| Scaffolds | 0 |
| Genes with eggNOG annotation | 18,777 (71.9%) |
| Genes with OrthoFinder ortholog | 13,826 (52.9%) |
| Genes with InterProScan annotation | 20,561 (78.7%) |
| Genes with RBH SwissProt hit | 7,077 (27.1%) |
| Genes with Tier 3 HMM domain hit | 6,213 (23.8%) |
| Genes with Tier 4 BLAST hit | 119 (0.5%) |
| Genes with >=1 TM helix (DeepTMHMM) | 4,739 (18.1%) |
| Genes with signal peptide SP (SignalP) | 2,992 (11.5%) |
| Genes with best human gene name | 16,938 (64.9%) |
| Genes with GO terms (any source) | 19,050 (72.9%) |
| Genes with KEGG pathway names | 7,848 (30.1%) |
| Median protein length | 378 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 6,641 | 25.4% |
| Silver | 7,581 | 29.0% |
| Bronze | 2,905 | 11.1% |
| Informative | 66 | 0.3% |
| Unclassified | 8,923 | 34.2% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 19,268 | 73.8% | Named human ortholog or informative eggNOG description |
| Partial | 3,166 | 12.1% | Named Pfam domain; no ortholog or functional description |
| Dark | 3,682 | 14.1% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 19,486 |
| Kinase | 875 |
| Transcription_Factor | 721 |
| Signalling | 620 |
| ECM_Adhesion | 571 |
| Innate_Immunity_PRR | 516 |
| GPCR | 497 |
| Innate_Immunity_Signalling | 428 |
| Symbiosis | 309 |
| Innate_Immunity_Effector | 269 |
| Ion_Channel | 257 |
| Neuron_Neuropeptide | 247 |
| Calcification_Universal | 233 |
| Lectin_Receptor | 185 |
| HSP_Chaperone | 182 |
| AMP_Venom | 160 |
| Innate_Immunity_Antiviral | 132 |
| Stem_Cell_Proliferation | 85 |
| Calcification_IonTransport | 82 |
| Stem_Cell_GMP | 69 |
| Calcification_Signalling | 64 |
| Tier4_SOMP | 51 |
| Tier4_SCRiP | 48 |
| Tier4_Neuropep | 19 |
| Stem_Cell_HSC | 10 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **27 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 18,176 | 69.6% |
| Receptor_candidate | 2,983 | 11.4% |
| Receptor | 3,072 | 11.8% |
| Ligand | 1,885 | 7.2% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **7,940 (30.4%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Xspp, human, mouse, *N. vectensis* |
| InterProScan | 5.78-109.0 | 17 member databases |
| DIAMOND | v2.2.1 | RBH vs SwissProt (Jan 2026) |
| SignalP | 6.0 | slow-sequential mode |
| DeepTMHMM | 1.0 | 2,000-sequence chunks |
| HMMER | 3.4 | 179 Tier 3 profiles |
| BioPython | v1.87 | pI and MW computation |
| 08_merge_annotate.py | | Integrates all 7 annotation layers; assigns confidence tiers, protein groups, CellChat roles, pI/MW, KEGG pathway names |

---

## Gene Lists

Pre-computed gene lists are in the `gene_lists/` directory:

- `Light_genes.txt` (19,268 genes)
- `Partial_genes.txt` (3,166 genes)
- `Dark_genes.txt` (3,682 genes)
- `CellChat_Ligand.txt` (1,885 genes)
- `CellChat_Receptor.txt` (3,072 genes)
- `CellChat_Receptor_candidate.txt` (2,983 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Xspp_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
