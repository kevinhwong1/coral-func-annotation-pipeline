# *Amil* Annotation Summary

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
| Genome version | Amil_scatlas_v1 |
| Proteome source | https://github.com/sebepedroslab/oculina-coral-sc-atlas |
| Proteome file | `Amil_long.pep.mod.fasta` |
| GFF3 file | `Amil_long.annot.gff3` |
| Pipeline run date | 2026-07-13 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Amil_master_annotation.tsv` (28,183 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 28,183 |
| Scaffolds | 588 |
| Genes with eggNOG annotation | 18,970 (67.3%) |
| Genes with OrthoFinder ortholog | 11,406 (40.5%) |
| Genes with InterProScan annotation | 18,074 (64.1%) |
| Genes with RBH SwissProt hit | 7,096 (25.2%) |
| Genes with Tier 3 HMM domain hit | 4,845 (17.2%) |
| Genes with Tier 4 BLAST hit | 52 (0.2%) |
| Genes with >=1 TM helix (DeepTMHMM) | 3,843 (13.6%) |
| Genes with signal peptide SP (SignalP) | 2,492 (8.8%) |
| Genes with best human gene name | 14,119 (50.1%) |
| Genes with GO terms (any source) | 15,925 (56.5%) |
| Genes with KEGG pathway names | 6,189 (22.0%) |
| Median protein length | 331 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 6,514 | 23.1% |
| Silver | 5,379 | 19.1% |
| Bronze | 2,232 | 7.9% |
| Informative | 24 | 0.1% |
| Unclassified | 14,034 | 49.8% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 18,098 | 64.2% | Named human ortholog or informative eggNOG description |
| Partial | 3,334 | 11.8% | Named Pfam domain; no ortholog or functional description |
| Dark | 6,751 | 24.0% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 23,084 |
| Kinase | 508 |
| ECM_Adhesion | 491 |
| Innate_Immunity_PRR | 485 |
| Signalling | 415 |
| Transcription_Factor | 398 |
| Innate_Immunity_Signalling | 376 |
| Lectin_Receptor | 266 |
| Ion_Channel | 262 |
| Calcification_Universal | 253 |
| Innate_Immunity_Effector | 247 |
| Symbiosis | 238 |
| GPCR | 208 |
| HSP_Chaperone | 186 |
| AMP_Venom | 179 |
| Neuron_Neuropeptide | 150 |
| Innate_Immunity_Antiviral | 140 |
| Stem_Cell_GMP | 85 |
| Stem_Cell_Proliferation | 60 |
| Calcification_IonTransport | 52 |
| Calcification_Signalling | 43 |
| Tier4_SCRiP | 35 |
| Tier4_Neuropep | 9 |
| Stem_Cell_HSC | 7 |
| Tier4_SOMP | 6 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **25 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 21,487 | 76.2% |
| Receptor_candidate | 2,669 | 9.5% |
| Receptor | 2,253 | 8.0% |
| Ligand | 1,774 | 6.3% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **6,696 (23.8%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Amil, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (18,098 genes)
- `Partial_genes.txt` (3,334 genes)
- `Dark_genes.txt` (6,751 genes)
- `CellChat_Ligand.txt` (1,774 genes)
- `CellChat_Receptor.txt` (2,253 genes)
- `CellChat_Receptor_candidate.txt` (2,669 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Amil_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
