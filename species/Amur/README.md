# *Amur* Annotation Summary

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
| Genome version | GCF_036669905.1 |
| Proteome source | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669905.1/ |
| Proteome file | `protein.faa` |
| GFF3 file | `genomic.gff` |
| Pipeline run date | 2026-07-20 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Amur_master_annotation.tsv` (42,312 genes x 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 42,312 |
| Scaffolds | 0 |
| Genes with eggNOG annotation | 30,281 (71.6%) |
| Genes with OrthoFinder ortholog | 16,651 (39.4%) |
| Genes with InterProScan annotation | 34,907 (82.5%) |
| Genes with RBH SwissProt hit | 8,009 (18.9%) |
| Genes with Tier 3 HMM domain hit | 12,126 (28.7%) |
| Genes with Tier 4 BLAST hit | 128 (0.3%) |
| Genes with >=1 TM helix (DeepTMHMM) | 10,556 (24.9%) |
| Genes with signal peptide SP (SignalP) | 6,948 (16.4%) |
| Genes with best human gene name | 22,292 (52.7%) |
| Genes with GO terms (any source) | 31,156 (73.6%) |
| Genes with KEGG pathway names | 11,651 (27.5%) |
| Median protein length | 435 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 7,317 | 17.3% |
| Silver | 9,912 | 23.4% |
| Bronze | 8,302 | 19.6% |
| Informative | 64 | 0.2% |
| Unclassified | 16,717 | 39.5% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 29,209 | 69.0% | Named human ortholog or informative eggNOG description |
| Partial | 7,591 | 17.9% | Named Pfam domain; no ortholog or functional description |
| Dark | 5,512 | 13.0% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 28,983 |
| GPCR | 2,569 |
| Innate_Immunity_PRR | 1,386 |
| ECM_Adhesion | 1,261 |
| Transcription_Factor | 1,059 |
| Kinase | 1,009 |
| Signalling | 876 |
| Innate_Immunity_Signalling | 682 |
| Ion_Channel | 662 |
| Calcification_Universal | 660 |
| Lectin_Receptor | 558 |
| Symbiosis | 429 |
| Innate_Immunity_Effector | 427 |
| Innate_Immunity_Antiviral | 310 |
| Neuron_Neuropeptide | 285 |
| AMP_Venom | 281 |
| HSP_Chaperone | 239 |
| Stem_Cell_GMP | 162 |
| Calcification_IonTransport | 124 |
| Stem_Cell_Proliferation | 117 |
| Calcification_Signalling | 97 |
| Tier4_SCRiP | 85 |
| Tier4_Neuropep | 33 |
| Stem_Cell_HSC | 11 |
| Tier4_SOMP | 7 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **126 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 25,621 | 60.6% |
| Receptor_candidate | 4,708 | 11.1% |
| Receptor | 7,974 | 18.8% |
| Ligand | 4,009 | 9.5% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **16,691 (39.4%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v3.1.0 | MSA mode, 4 species: Amur, human, mouse, *N. vectensis* |
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

- `Light_genes.txt` (29,209 genes)
- `Partial_genes.txt` (7,591 genes)
- `Dark_genes.txt` (5,512 genes)
- `CellChat_Ligand.txt` (4,009 genes)
- `CellChat_Receptor.txt` (7,974 genes)
- `CellChat_Receptor_candidate.txt` (4,708 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Amur_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
