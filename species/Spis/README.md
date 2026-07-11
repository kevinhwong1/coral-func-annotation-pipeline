# *Spis* Annotation Summary

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
| Genome version | Spis_scatlas_v1_Levy2021 |
| Proteome source | https://github.com/sebepedroslab/oculina-coral-sc-atlas |
| Proteome file | `Spis_long.pep.mod.fasta` |
| GFF3 file | `Spis_long.annot.gff3` |
| Pipeline run date | 2026-07-11 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Spis_master_annotation.tsv` (29,419 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 29,419 |
| Scaffolds | 3,122 |
| Genes with eggNOG annotation | 19,286 (65.6%) |
| Genes with OrthoFinder ortholog | 11,998 (40.8%) |
| Genes with InterProScan annotation | 22,581 (76.8%) |
| Genes with RBH SwissProt hit | 8,553 (29.1%) |
| Genes with Tier 3 HMM domain hit | 7,060 (24.0%) |
| Genes with Tier 4 BLAST hit | 102 (0.3%) |
| Genes with >=1 TM helix (DeepTMHMM) | 5,908 (20.1%) |
| Genes with signal peptide SP (SignalP) | 3,622 (12.3%) |
| Genes with best human gene name | 15,832 (53.8%) |
| Genes with GO terms (any source) | 19,859 (67.5%) |
| Genes with KEGG pathway names | 7,023 (23.9%) |
| Median protein length | 362 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 7,416 | 25.2% |
| Silver | 5,564 | 18.9% |
| Bronze | 4,288 | 14.6% |
| Informative | 53 | 0.2% |
| Unclassified | 12,098 | 41.1% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 20,073 | 68.2% | Named human ortholog or informative eggNOG description |
| Partial | 4,564 | 15.5% | Named Pfam domain; no ortholog or functional description |
| Dark | 4,782 | 16.3% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 21,791 |
| GPCR | 1,297 |
| ECM_Adhesion | 744 |
| Innate_Immunity_PRR | 614 |
| Kinase | 613 |
| Transcription_Factor | 528 |
| Signalling | 506 |
| Calcification_Universal | 468 |
| Innate_Immunity_Signalling | 423 |
| Lectin_Receptor | 357 |
| Innate_Immunity_Effector | 341 |
| Ion_Channel | 333 |
| Symbiosis | 279 |
| Innate_Immunity_Antiviral | 244 |
| HSP_Chaperone | 179 |
| AMP_Venom | 175 |
| Neuron_Neuropeptide | 157 |
| Stem_Cell_GMP | 76 |
| Tier4_SCRiP | 75 |
| Stem_Cell_Proliferation | 71 |
| Calcification_IonTransport | 62 |
| Calcification_Signalling | 53 |
| Tier4_SOMP | 12 |
| Tier4_Neuropep | 11 |
| Stem_Cell_HSC | 10 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **79 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 19,859 | 67.5% |
| Receptor_candidate | 3,153 | 10.7% |
| Receptor | 4,028 | 13.7% |
| Ligand | 2,379 | 8.1% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **9,560 (32.5%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Spis, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (20,073 genes)
- `Partial_genes.txt` (4,564 genes)
- `Dark_genes.txt` (4,782 genes)
- `CellChat_Ligand.txt` (2,379 genes)
- `CellChat_Receptor.txt` (4,028 genes)
- `CellChat_Receptor_candidate.txt` (3,153 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Spis_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
