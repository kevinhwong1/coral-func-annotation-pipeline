# *Galaxea fascicularis* (Gfas) Annotation Summary

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
| Genome version | Gfas_v1.0 |
| Proteome source | [gfas.reefgenomics.org](http://gfas.reefgenomics.org/) |
| Proteome file | `gfas_1.0.proteins.fasta` |
| GFF3 file | `gfas_1.0.genes.gff3` |
| Pipeline run date | 2026-06-05 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Gfas_master_annotation.tsv` (22,418 genes × 52 columns; `protein_sequence` excluded) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.


---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 22,418 |
| Scaffolds | 6,089 |
| Genes with eggNOG annotation | 15,381 (68.6%) |
| Genes with OrthoFinder ortholog | 17,803 (79.4%) |
| Genes with InterProScan annotation | 17,590 (78.5%) |
| Genes with RBH SwissProt hit | 7,618 (34.0%) |
| Genes with Tier 3 HMM domain hit | 5,224 (23.3%) |
| Genes with Tier 4 BLAST hit | 57 (0.3%) |
| Genes with ≥1 TM helix (DeepTMHMM) | 4,168 (18.6%) |
| Genes with signal peptide SP (SignalP) | 2,355 (10.5%) |
| Genes with best human gene name | 14,074 (62.8%) |
| Genes with GO terms (any source) | 15,147 (67.6%) |
| Genes with KEGG pathway names | 6,049 (27.0%) |
| Median protein length | 314 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 6,996 | 31.2% |
| Silver | 5,153 | 23.0% |
| Bronze | 2,398 | 10.7% |
| Informative | 27 | 0.1% |
| Unclassified | 7,844 | 35.0% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 16,060 | 71.6% | Named human ortholog or informative eggNOG description |
| Partial | 2,230 | 9.9% | Named Pfam domain; no ortholog or functional description |
| Dark | 4,128 | 18.4% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 16,863 |
| GPCR | 789 |
| ECM_Adhesion | 506 |
| Transcription_Factor | 430 |
| Kinase | 420 |
| Signalling | 407 |
| Innate_Immunity_PRR | 333 |
| Calcification_Universal | 328 |
| Innate_Immunity_Signalling | 313 |
| Lectin_Receptor | 307 |
| Ion_Channel | 266 |
| Innate_Immunity_Antiviral | 247 |
| Innate_Immunity_Effector | 246 |
| Symbiosis | 232 |
| HSP_Chaperone | 158 |
| Neuron_Neuropeptide | 150 |
| AMP_Venom | 144 |
| Stem_Cell_GMP | 64 |
| Stem_Cell_Proliferation | 62 |
| Calcification_IonTransport | 57 |
| Calcification_Signalling | 35 |
| Tier4_SCRiP | 31 |
| Tier4_Neuropep | 12 |
| Tier4_SOMP | 11 |
| Stem_Cell_HSC | 7 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **48 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 15,656 | 69.8% |
| Receptor_candidate | 2,551 | 11.4% |
| Receptor | 2,530 | 11.3% |
| Ligand | 1,681 | 7.5% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **6,762 (30.2%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v2 | MSA mode, 4 species: Gfas, human, mouse, *N. vectensis* |
| InterProScan | 5.78-109.0 | 17 member databases |
| DIAMOND | v2.2.1 | RBH vs SwissProt (Jan 2026) |
| SignalP | 6.0 | slow-sequential mode |
| DeepTMHMM | 1.0 | 14 chunks × 2,000 sequences |
| HMMER | 3.4 | 179 Tier 3 profiles |
| BioPython | v1.87 | pI and MW computation |
| 08_merge_annotate.py | — | Integrates all 7 annotation layers; assigns confidence tiers, protein groups, CellChat roles, pI/MW, KEGG pathway names |

---

## Gene Lists

Pre-computed gene lists are in the `gene_lists/` directory:

- `Light_genes.txt` (16,060 genes)
- `Partial_genes.txt` (2,230 genes)
- `Dark_genes.txt` (4,128 genes)
- `CellChat_Ligand.txt` (1,681 genes)
- `CellChat_Receptor.txt` (2,530 genes)
- `CellChat_Receptor_candidate.txt` (2,551 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Gfas_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
