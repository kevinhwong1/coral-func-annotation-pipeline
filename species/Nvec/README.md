# *Nvec* Annotation Summary

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
| Genome version | NV2g.20240221 |
| Proteome source | https://simrbase.stowers.org/nematostella |
| Proteome file | `NV2g.20240221.protein.fa` |
| GFF3 file | `NV2g.20240221.gff` |
| Pipeline run date | 2026-07-06 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Nvec_master_annotation.tsv` (32,495 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 32,495 |
| Scaffolds | 29 |
| Genes with eggNOG annotation | 21,106 (65.0%) |
| Genes with OrthoFinder ortholog | 14,170 (43.6%) |
| Genes with InterProScan annotation | 23,091 (71.1%) |
| Genes with RBH SwissProt hit | 7,838 (24.1%) |
| Genes with Tier 3 HMM domain hit | 7,359 (22.6%) |
| Genes with Tier 4 BLAST hit | 95 (0.3%) |
| Genes with >=1 TM helix (DeepTMHMM) | 6,122 (18.8%) |
| Genes with signal peptide SP (SignalP) | 4,158 (12.8%) |
| Genes with best human gene name | 17,816 (54.8%) |
| Genes with GO terms (any source) | 21,921 (67.5%) |
| Genes with KEGG pathway names | 8,846 (27.2%) |
| Median protein length | 339 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 6,992 | 21.5% |
| Silver | 8,006 | 24.6% |
| Bronze | 3,888 | 12.0% |
| Informative | 43 | 0.1% |
| Unclassified | 13,566 | 41.7% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 21,344 | 65.7% | Named human ortholog or informative eggNOG description |
| Partial | 3,618 | 11.1% | Named Pfam domain; no ortholog or functional description |
| Dark | 7,533 | 23.2% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 24,805 |
| GPCR | 854 |
| Transcription_Factor | 737 |
| ECM_Adhesion | 705 |
| Signalling | 627 |
| Kinase | 588 |
| Innate_Immunity_PRR | 475 |
| Innate_Immunity_Signalling | 474 |
| Calcification_Universal | 385 |
| Lectin_Receptor | 384 |
| Symbiosis | 381 |
| Ion_Channel | 380 |
| Neuron_Neuropeptide | 318 |
| Innate_Immunity_Effector | 312 |
| AMP_Venom | 268 |
| HSP_Chaperone | 174 |
| Innate_Immunity_Antiviral | 156 |
| Stem_Cell_Proliferation | 130 |
| Calcification_IonTransport | 96 |
| Stem_Cell_GMP | 88 |
| Tier4_SCRiP | 60 |
| Calcification_Signalling | 46 |
| Tier4_Neuropep | 26 |
| Stem_Cell_HSC | 20 |
| Tier4_SOMP | 6 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **57 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 22,292 | 68.6% |
| Receptor_candidate | 3,806 | 11.7% |
| Receptor | 3,570 | 11.0% |
| Ligand | 2,827 | 8.7% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **10,203 (31.4%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 3 species: Nvec, human, mouse |
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

- `Light_genes.txt` (21,344 genes)
- `Partial_genes.txt` (3,618 genes)
- `Dark_genes.txt` (7,533 genes)
- `CellChat_Ligand.txt` (2,827 genes)
- `CellChat_Receptor.txt` (3,570 genes)
- `CellChat_Receptor_candidate.txt` (3,806 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Nvec_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
