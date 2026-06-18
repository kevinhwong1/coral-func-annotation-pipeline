# *Pocillopora damicornis* (Pdam) Annotation Summary

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
| Genome version | Pdam_v1.0 |
| Proteome source | [pdam.reefgenomics.org](http://pdam.reefgenomics.org/) |
| Proteome file | `pdam_proteins.fasta` |
| GFF3 file | `pdam_genome.gff3` |
| Pipeline run date | 2026-06-16 |
| Merge script | `08_merge_annotate.py` |
| Master table | `Pdam_master_annotation.tsv` (26,077 genes × 53 columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | 26,077 |
| Scaffolds | 1,968 |
| Genes with eggNOG annotation | 15,370 (58.9%) |
| Genes with OrthoFinder ortholog | 19,908 (76.3%) |
| Genes with InterProScan annotation | 17,943 (68.8%) |
| Genes with RBH SwissProt hit | 7,659 (29.4%) |
| Genes with Tier 3 HMM domain hit | 5,381 (20.6%) |
| Genes with Tier 4 BLAST hit | 68 (0.3%) |
| Genes with ≥1 TM helix (DeepTMHMM) | 5,157 (19.8%) |
| Genes with signal peptide SP (SignalP) | 2,502 (9.6%) |
| Genes with best human gene name | 14,481 (55.5%) |
| Genes with GO terms (any source) | 15,486 (59.4%) |
| Genes with KEGG pathway names | 6,080 (23.3%) |
| Median protein length | 306 aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | 7,066 | 27.1% |
| Silver | 5,620 | 21.6% |
| Bronze | 2,630 | 10.1% |
| Informative | 15 | 0.1% |
| Unclassified | 10,746 | 41.2% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | 16,756 | 64.3% | Named human ortholog or informative eggNOG description |
| Partial | 2,435 | 9.3% | Named Pfam domain; no ortholog or functional description |
| Dark | 6,886 | 26.4% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
| Unclassified | 20,134 |
| GPCR | 1,074 |
| ECM_Adhesion | 550 |
| Transcription_Factor | 465 |
| Kinase | 419 |
| Signalling | 408 |
| Innate_Immunity_PRR | 363 |
| Innate_Immunity_Signalling | 353 |
| Calcification_Universal | 326 |
| Lectin_Receptor | 289 |
| Ion_Channel | 283 |
| Symbiosis | 254 |
| Innate_Immunity_Effector | 249 |
| AMP_Venom | 176 |
| Innate_Immunity_Antiviral | 155 |
| Neuron_Neuropeptide | 145 |
| HSP_Chaperone | 132 |
| Stem_Cell_Proliferation | 66 |
| Stem_Cell_GMP | 58 |
| Calcification_IonTransport | 56 |
| Calcification_Signalling | 51 |
| Tier4_SCRiP | 45 |
| Tier4_SOMP | 10 |
| Tier4_Neuropep | 10 |
| Stem_Cell_HSC | 6 |

Dual-role genes with a secondary group assignment (`protein_group_secondary`): **57 genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | 18,300 | 70.2% |
| Receptor_candidate | 3,086 | 11.8% |
| Receptor | 2,952 | 11.3% |
| Ligand | 1,739 | 6.7% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **7,777 (29.8%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v2 | MSA mode, 4 species: Pdam, human, mouse, *N. vectensis* |
| InterProScan | 5.78-109.0 | 17 member databases |
| DIAMOND | v2.2.1 | RBH vs SwissProt (Jan 2026) |
| SignalP | 6.0 | slow-sequential mode |
| DeepTMHMM | 1.0 | 2,000-sequence chunks |
| HMMER | 3.4 | 179 Tier 3 profiles |
| BioPython | v1.87 | pI and MW computation |
| 08_merge_annotate.py | — | Integrates all 7 annotation layers; assigns confidence tiers, protein groups, CellChat roles, pI/MW, KEGG pathway names |

---

## Gene Lists

Pre-computed gene lists are in the `gene_lists/` directory:

- `Light_genes.txt` (16,756 genes)
- `Partial_genes.txt` (2,435 genes)
- `Dark_genes.txt` (6,886 genes)
- `CellChat_Ligand.txt` (1,739 genes)
- `CellChat_Receptor.txt` (2,952 genes)
- `CellChat_Receptor_candidate.txt` (3,086 genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`Pdam_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
