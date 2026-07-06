# *Ofav* Annotation Summary

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
| Genome version | Ofav_gen_17 |
| Proteome source | University of Miami long-read assembly (internal) |
| Proteome file | `Orbicella_faveolata_gen_17.proteins.fa` |
| GFF3 file | `Orbicella_faveolata_gen_17.gff3` |
| Pipeline run date | 2026-06-24 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Ofav_master_annotation.tsv` (32,172 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 32,172 |
| Scaffolds | 35 |
| Genes with eggNOG annotation | 19,606 (60.9%) |
| Genes with OrthoFinder ortholog | 12,852 (39.9%) |
| Genes with InterProScan annotation | 22,076 (68.6%) |
| Genes with RBH SwissProt hit | 8,085 (25.1%) |
| Genes with Tier 3 HMM domain hit | 7,123 (22.1%) |
| Genes with Tier 4 BLAST hit | 89 (0.3%) |
| Genes with >=1 TM helix (DeepTMHMM) | 6,864 (21.3%) |
| Genes with signal peptide SP (SignalP) | 4,382 (13.6%) |
| Genes with best human gene name | 16,705 (51.9%) |
| Genes with GO terms (any source) | 20,179 (62.7%) |
| Genes with KEGG pathway names | 8,034 (25.0%) |
| Median protein length | 335 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 7,370 | 22.9% |
| Silver | 6,052 | 18.8% |
| Bronze | 4,314 | 13.4% |
| Informative | 36 | 0.1% |
| Unclassified | 14,400 | 44.8% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 19,972 | 62.1% | Named human ortholog or informative eggNOG description |
| Partial | 4,150 | 12.9% | Named Pfam domain; no ortholog or functional description |
| Dark | 8,050 | 25.0% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 24,386 |
| GPCR | 1,371 |
| ECM_Adhesion | 674 |
| Transcription_Factor | 608 |
| Kinase | 590 |
| Calcification_Universal | 582 |
| Signalling | 534 |
| Innate_Immunity_PRR | 474 |
| Innate_Immunity_Signalling | 447 |
| Ion_Channel | 403 |
| Innate_Immunity_Effector | 356 |
| Symbiosis | 325 |
| Lectin_Receptor | 308 |
| AMP_Venom | 230 |
| HSP_Chaperone | 185 |
| Neuron_Neuropeptide | 179 |
| Innate_Immunity_Antiviral | 160 |
| Stem_Cell_Proliferation | 75 |
| Stem_Cell_GMP | 71 |
| Calcification_IonTransport | 63 |
| Tier4_SCRiP | 58 |
| Calcification_Signalling | 55 |
| Tier4_SOMP | 15 |
| Tier4_Neuropep | 13 |
| Stem_Cell_HSC | 10 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **74 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 21,081 | 65.5% |
| Receptor_candidate | 3,969 | 12.3% |
| Receptor | 4,014 | 12.5% |
| Ligand | 3,108 | 9.7% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **11,091 (34.5%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Ofav, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (19,972 genes)
- `Partial_genes.txt` (4,150 genes)
- `Dark_genes.txt` (8,050 genes)
- `CellChat_Ligand.txt` (3,108 genes)
- `CellChat_Receptor.txt` (4,014 genes)
- `CellChat_Receptor_candidate.txt` (3,969 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Ofav_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
