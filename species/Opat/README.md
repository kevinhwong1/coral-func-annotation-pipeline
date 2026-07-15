# *Opat* Annotation Summary

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
| Genome version | Ocupat_scatlas_v1 |
| Proteome source | https://github.com/sebepedroslab/oculina-coral-sc-atlas |
| Proteome file | `Ocupat_long.pep.fasta` |
| GFF3 file | `Ocupat_long.annot.gff3` |
| Pipeline run date | 2026-07-14 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Ocupat_master_annotation.tsv` (39,482 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 39,482 |
| Scaffolds | 14 |
| Genes with eggNOG annotation | 24,245 (61.4%) |
| Genes with OrthoFinder ortholog | 15,944 (40.4%) |
| Genes with InterProScan annotation | 26,899 (68.1%) |
| Genes with RBH SwissProt hit | 7,981 (20.2%) |
| Genes with Tier 3 HMM domain hit | 8,812 (22.3%) |
| Genes with Tier 4 BLAST hit | 137 (0.3%) |
| Genes with >=1 TM helix (DeepTMHMM) | 7,666 (19.4%) |
| Genes with signal peptide SP (SignalP) | 5,853 (14.8%) |
| Genes with best human gene name | 20,215 (51.2%) |
| Genes with GO terms (any source) | 23,341 (59.1%) |
| Genes with KEGG pathway names | 8,781 (22.2%) |
| Median protein length | 367 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 7,260 | 18.4% |
| Silver | 9,215 | 23.3% |
| Bronze | 5,442 | 13.8% |
| Informative | 87 | 0.2% |
| Unclassified | 17,478 | 44.3% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 25,444 | 64.4% | Named human ortholog or informative eggNOG description |
| Partial | 5,058 | 12.8% | Named Pfam domain; no ortholog or functional description |
| Dark | 8,980 | 22.7% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 29,851 |
| GPCR | 1,210 |
| ECM_Adhesion | 1,179 |
| Calcification_Universal | 1,002 |
| Transcription_Factor | 893 |
| Innate_Immunity_PRR | 801 |
| Kinase | 563 |
| Innate_Immunity_Signalling | 557 |
| Lectin_Receptor | 554 |
| Signalling | 500 |
| Innate_Immunity_Effector | 435 |
| Ion_Channel | 376 |
| Symbiosis | 343 |
| AMP_Venom | 285 |
| Innate_Immunity_Antiviral | 212 |
| Neuron_Neuropeptide | 166 |
| HSP_Chaperone | 152 |
| Tier4_SCRiP | 91 |
| Stem_Cell_Proliferation | 83 |
| Stem_Cell_GMP | 64 |
| Calcification_Signalling | 59 |
| Calcification_IonTransport | 57 |
| Tier4_SOMP | 28 |
| Tier4_Neuropep | 13 |
| Stem_Cell_HSC | 8 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **95 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 26,122 | 66.2% |
| Receptor_candidate | 4,625 | 11.7% |
| Receptor | 4,515 | 11.4% |
| Ligand | 4,220 | 10.7% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **13,360 (33.8%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Ocupat, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (25,444 genes)
- `Partial_genes.txt` (5,058 genes)
- `Dark_genes.txt` (8,980 genes)
- `CellChat_Ligand.txt` (4,220 genes)
- `CellChat_Receptor.txt` (4,515 genes)
- `CellChat_Receptor_candidate.txt` (4,625 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Ocupat_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
