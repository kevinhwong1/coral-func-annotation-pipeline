# *Acer* Annotation Summary

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
| Genome version | GCA_032359415.1 |
| Proteome source | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/ |
| Proteome file | `acer_protein.faa` |
| GFF3 file | `genomic.gff` |
| Pipeline run date | 2026-06-18 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Acer_master_annotation.tsv` (28,059 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 28,059 |
| Scaffolds | 379 |
| Genes with eggNOG annotation | 17,294 (61.6%) |
| Genes with OrthoFinder ortholog | 11,486 (40.9%) |
| Genes with InterProScan annotation | 16,859 (60.1%) |
| Genes with RBH SwissProt hit | 7,486 (26.7%) |
| Genes with Tier 3 HMM domain hit | 4,623 (16.5%) |
| Genes with Tier 4 BLAST hit | 48 (0.2%) |
| Genes with >=1 TM helix (DeepTMHMM) | 4,251 (15.2%) |
| Genes with signal peptide SP (SignalP) | 2,193 (7.8%) |
| Genes with best human gene name | 14,144 (50.4%) |
| Genes with GO terms (any source) | 15,259 (54.4%) |
| Genes with KEGG pathway names | 5,839 (20.8%) |
| Median protein length | 295 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 6,846 | 24.4% |
| Silver | 5,186 | 18.5% |
| Bronze | 2,339 | 8.3% |
| Informative | 14 | 0.0% |
| Unclassified | 13,674 | 48.7% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 17,392 | 62.0% | Named human ortholog or informative eggNOG description |
| Partial | 2,807 | 10.0% | Named Pfam domain; no ortholog or functional description |
| Dark | 7,860 | 28.0% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 23,071 |
| GPCR | 756 |
| Transcription_Factor | 469 |
| Signalling | 410 |
| ECM_Adhesion | 377 |
| Kinase | 369 |
| Innate_Immunity_Signalling | 332 |
| Calcification_Universal | 255 |
| Innate_Immunity_PRR | 255 |
| Ion_Channel | 251 |
| Symbiosis | 244 |
| Innate_Immunity_Effector | 214 |
| Lectin_Receptor | 196 |
| Innate_Immunity_Antiviral | 168 |
| HSP_Chaperone | 149 |
| Neuron_Neuropeptide | 148 |
| AMP_Venom | 123 |
| Stem_Cell_GMP | 76 |
| Stem_Cell_Proliferation | 57 |
| Calcification_IonTransport | 47 |
| Calcification_Signalling | 41 |
| Tier4_SCRiP | 31 |
| Tier4_Neuropep | 9 |
| Stem_Cell_HSC | 6 |
| Tier4_SOMP | 5 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **51 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 21,583 | 76.9% |
| Receptor_candidate | 2,658 | 9.5% |
| Receptor | 2,278 | 8.1% |
| Ligand | 1,540 | 5.5% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **6,476 (23.1%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v2 | MSA mode, 4 species: Acer, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (17,392 genes)
- `Partial_genes.txt` (2,807 genes)
- `Dark_genes.txt` (7,860 genes)
- `CellChat_Ligand.txt` (1,540 genes)
- `CellChat_Receptor.txt` (2,278 genes)
- `CellChat_Receptor_candidate.txt` (2,658 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Acer_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
