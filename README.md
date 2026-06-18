# coral-func-annotation-pipeline

A pipeline for functional annotation of coral proteomes. Each gene model is annotated with orthology assignments, domain architecture, signal peptide predictions, transmembrane topology, and functional classification. Outputs a single master annotation table per species for use in comparative genomics, single-cell RNA-seq, and cell communication analyses.

Annotation status tiers (Light/Partial/Dark) follow the framework of [Stephens et al. 2026](#references).

---

## Contents

- [Species](#species)
- [Pipeline Overview](#pipeline-overview)
- [Output](#output)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Classification Summary](#classification-summary)
- [Tool Versions](#tool-versions)
- [Parameters](#parameters)
- [Abbreviations](#abbreviations)
- [References](#references)
- [Citation](#citation)
- [License](#license)

---

## Species

### Genome Sources

Proteome sequences are not hosted in this repository. Download from the sources below before running the pipeline. Proteomes marked with * must be downloaded manually and transferred to the cluster (direct wget is not available from the source).

| Species | Abbrev | Source | Proteome file |
|---------|--------|--------|---------------|
| *Galaxea fascicularis* | Gfas | [gfas.reefgenomics.org](http://gfas.reefgenomics.org/) | `gfas_1.0.proteins.fasta` |
| *Pocillopora damicornis* | Pdam | [pdam.reefgenomics.org](http://pdam.reefgenomics.org/) | `pdam_proteins.fasta` |
| *Acropora millepora* | Amil | [sebepedroslab/oculina-coral-sc-atlas](https://github.com/sebepedroslab/oculina-coral-sc-atlas/blob/master/data/reference/Amil_long.pep.fasta) * | `Amil_long.pep.fasta` |
| *Stylophora pistillata* | Spis | [sebepedroslab/oculina-coral-sc-atlas](https://github.com/sebepedroslab/oculina-coral-sc-atlas/blob/master/data/reference/Spis_long.pep.fasta) * | `Spis_long.pep.fasta` |
| *Oculina patagonica* | Opat | [sebepedroslab/oculina-coral-sc-atlas](https://github.com/sebepedroslab/oculina-coral-sc-atlas/blob/master/data/reference/Ocupat_long.pep.fasta) * | `Ocupat_long.pep.fasta` |
| *Acropora cervicornis* | Acer | [NCBI GCA_032359415.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/) | — |
| *Orbicella faveolata* | Ofav | [NCBI GCA_042242905.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_042242905.1/) | — |
| *Xenia* spp. | Xspp | [Carnegie Endosymbiosis](https://cmo.carnegiescience.edu/endosymbiosis/genome/) | `xenSp1.proteins.fa` |
| *Acropora muricata* | Amur | [NCBI GCF_036669905.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669905.1/) | — |
| *Nematostella vectensis* | Nvec | [SimrBase](https://simrbase.stowers.org/nematostella) | `NV2g.20240221.protein.fa` |

> For NCBI genomes, download the predicted protein FASTA (`*_protein.faa`) using the NCBI Datasets CLI: `datasets download genome accession <GCA_ID> --include protein`. Proteomes from GitHub (Amil, Spis, Opat) must be downloaded via the GitHub interface or `git clone` and transferred to the cluster manually.
>
> *Nematostella vectensis* (Nvec) is used as an outgroup in OrthoFinder and will be annotated in a future pipeline run.

### Annotation Results

| Species | Abbrev | n genes | Gold | Silver | Bronze | Informative | Unclassified | Light | Partial | Dark | Run date |
|---------|--------|---------|------|--------|--------|-------------|--------------|-------|---------|------|----------|
| *Galaxea fascicularis* | Gfas | 22,418 | 31.2% | 23.0% | 10.7% | 0.1% | 35.0% | 71.6% | 9.9% | 18.4% | 2026-06-05 |
| *Pocillopora damicornis* | Pdam | 26,077 | | | | | | | | | in progress |
| *Acropora millepora* | Amil | | | | | | | | | | planned |
| *Stylophora pistillata* | Spis | | | | | | | | | | planned |
| *Oculina patagonica* | Opat | | | | | | | | | | planned |
| *Acropora cervicornis* | Acer | | | | | | | | | | planned |
| *Orbicella faveolata* | Ofav | | | | | | | | | | planned |
| *Xenia* spp. | Xspp | | | | | | | | | | planned |
| *Acropora muricata* | Amur | | | | | | | | | | planned |
| *Nematostella vectensis* | Nvec | | | | | | | | | | planned |

---

## Pipeline Overview

Eight jobs are submitted per species via LSF/BSUB. BSUB job scripts are generated at runtime by `submit_pipeline.sh` and written to the species run directory. They are not stored as standalone files.

| Job | Tool | Purpose | Queue | Wall time |
|-----|------|---------|-------|-----------|
| 01 | eggNOG-mapper v2 | Orthology, GO, KEGG, COG | bigmem | 12 hr |
| 02 | OrthoFinder v2 | Phylogenetic orthologs (human, mouse, *N. vectensis*) | bigmem | 120 hr |
| 03 | InterProScan 5.78-109.0 | Domain annotation (17 databases) | bigmem | 48 hr |
| 04 | DIAMOND v2.2.1 | Reciprocal best hits vs SwissProt | bigmem | 8 hr |
| 05 | SignalP 6.0 | Signal peptide prediction | general | 20 hr |
| 05b | DeepTMHMM 1.0 | Transmembrane topology (2,000-seq chunks) | general | 12 hr/chunk |
| 06 | HMMER 3.4 | Tier 3 custom HMM profiles (179 profiles) | general | 6 hr |
| 07 | DIAMOND v2.2.1 | Tier 4 curated coral BLAST (17 sequences) | general | 1 hr |
| 08 | 08_merge_annotate.py | Merge all layers → master annotation table | general | 2 hr |

See [`pipeline/README.md`](pipeline/README.md) for full submission instructions, configuration, and HPC notes.

---

## Output

Each species produces a **master annotation table** (53 columns per gene, including protein sequences; see [`docs/column_descriptions.md`](docs/column_descriptions.md) for the full schema) and a set of gene lists. The `protein_sequence` column (column 40) contains the full amino acid sequence for each gene model.

See [`docs/classification_system.md`](docs/classification_system.md) for full classification logic.

---

## Repository Structure

```
coral-func-annotation-pipeline/
│
├── README.md                            # This file
├── LICENSE
│
├── pipeline/
│   ├── README.md                        # Submission instructions, HPC notes, job config
│   ├── 00_setup_annotation_env.sh       # One-time environment and database setup
│   ├── submit_pipeline.sh               # Master job submission script (generates BSUB scripts)
│   ├── run_config.sh                    # Species paths and pipeline variables
│   └── 08_merge_annotate.py             # Merge and classification script (~1,700 lines)
│
├── databases/
│   └── README.md                        # Database names, versions, download sources
│
├── species/
│   ├── Gfas/                            # Complete
│   │   ├── Gfas_annotation_summary.md       # Results summary (tiers, groups, CellChat)
│   │   ├── Gfas_master_annotation.tsv        # 53-column table (includes protein_sequence)
│   │   │                                     # See docs/column_descriptions.md for schema
│   │   ├── intermediate/
│   │   │   ├── eggnog_results.tsv
│   │   │   ├── interproscan_results.tsv
│   │   │   ├── rbh_swissprot.tsv
│   │   │   ├── signalp_summary.tsv
│   │   │   ├── tier3_hmm_hits.tsv
│   │   │   ├── tier4_blast_hits.tsv
│   │   │   ├── deeptmhmm_merged.gff3
│   │   │   └── orthofinder/
│   │   │       ├── Orthogroups.tsv
│   │   │       ├── Statistics_Overall.tsv
│   │   │       └── Orthologues/
│   │   └── gene_lists/
│   │       ├── Light_genes.txt
│   │       ├── Partial_genes.txt
│   │       ├── Dark_genes.txt
│   │       ├── CellChat_Ligand.txt
│   │       ├── CellChat_Receptor.txt
│   │       ├── CellChat_Receptor_candidate.txt
│   │       └── by_protein_group/        # One .txt file per functional group
│   ├── Pdam/                            # Same structure as Gfas
│   ├── Amil/                            # Same structure as Gfas
│   ├── Spis/                            # Same structure as Gfas
│   ├── Opat/                            # Same structure as Gfas
│   ├── Acer/                            # Same structure as Gfas
│   ├── Ofav/                            # Same structure as Gfas
│   ├── Xspp/                            # Same structure as Gfas
│   ├── Amur/                            # Same structure as Gfas
│   └── Nvec/                            # Same structure as Gfas
│
├── docs/
│   ├── column_descriptions.md           # Full 53-column schema with types, examples, citations
│   └── classification_system.md         # Tier logic, protein groups, CellChat roles, citations
│
└── examples/
    ├── load_annotation.R                # Quick Seurat/tidyverse join example
    └── load_annotation.py               # Quick pandas example
```

---

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/kevinhwong1/coral-func-annotation-pipeline.git
cd coral-func-annotation-pipeline

# 2. One-time environment and database setup
# See pipeline/README.md for prerequisites (SignalP license, DeepTMHMM academic license)
bash pipeline/00_setup_annotation_env.sh

# 3. Configure your species run
# Edit the paths at the top of run_config.sh for your species proteome and GFF3
cp pipeline/run_config.sh my_run_config.sh

# 4. Submit pipeline (LSF/BSUB)
bash pipeline/submit_pipeline.sh <SPECIES>
# e.g. bash pipeline/submit_pipeline.sh Pdam
```

See [`pipeline/README.md`](pipeline/README.md) for detailed instructions including database setup, DeepTMHMM chunking, known issues, and HPC queue configuration.

---

## Classification Summary

Genes are classified along four independent axes. Full logic is in [`docs/classification_system.md`](docs/classification_system.md).

**Confidence tier**: evidence strength from orthology and homology sources.
`Gold` > `Silver` > `Bronze` > `Informative` > `Unclassified`

**Annotation status**: functional characterisation in public databases, following Stephens et al. 2026.
`Light` (named ortholog or description) | `Partial` (domain only) | `Dark` (no evidence)

**Protein group**: 24 functional groups for use in Seurat and CellChat analyses.

**CellChat role**: `Ligand`, `Receptor`, `Receptor_candidate`, or `Neither`.

---

## Tool Versions

| Tool | Version | Publication |
|------|---------|-------------|
| eggNOG-mapper | v2 | [Cantalapiedra et al. 2021](https://doi.org/10.1093/molbev/msab293) |
| OrthoFinder | v2 | [Emms & Kelly 2019](https://doi.org/10.1186/s13059-019-1832-y) |
| InterProScan | 5.78-109.0 | [Jones et al. 2014](https://doi.org/10.1093/bioinformatics/btu031) |
| DIAMOND | v2.2.1 | [Buchfink et al. 2021](https://doi.org/10.1038/s41592-021-01101-x) |
| SignalP | 6.0 | [Teufel et al. 2022](https://doi.org/10.1038/s41587-021-01156-3) |
| DeepTMHMM | 1.0 | [Hallgren et al. 2022](https://doi.org/10.1101/2022.04.08.487609) |
| HMMER | 3.4 | [Eddy 2011](https://doi.org/10.1371/journal.pcbi.1002195) |
| BioPython | v1.87 | [Cock et al. 2009](https://doi.org/10.1093/bioinformatics/btp163) |

---

## Parameters

Parameters for each annotation step. Tools were run with default parameters unless otherwise noted.

| Step | Tool | Parameters |
|------|------|------------|
| eggNOG orthology | eggNOG-mapper v2 | Default sensitivity (`--sensmode default`); eggNOG 5.0 database |
| Ortholog inference | OrthoFinder v2 | MSA mode (`-M msa`); 4 species: coral, human, mouse, *N. vectensis* |
| Domain annotation | InterProScan 5.78-109.0 | Per-database cutoffs as distributed (Pfam uses GA bit-score thresholds; PANTHER e-value ≤ 1e-3; other databases use database-specific defaults). No user-level e-value filter applied. |
| RBH vs SwissProt | BLAST+ blastp | e-value ≤ 1e-5; `-max_target_seqs 1` (applied in both directions independently) |
| Signal peptide | SignalP 6.0 | `--mode slow-sequential`; `--organism eukarya`; classification threshold set by model output (no user-adjustable cutoff) |
| TM topology | DeepTMHMM 1.0 | Default model parameters; proteomes split into 2,000-sequence chunks; CPU inference |
| Tier 3 HMM | HMMER 3.4 hmmsearch | `-E 1e-3 --domE 1e-3`; 179 custom coral/metazoan HMM profiles |
| Tier 4 BLAST | BLAST+ blastp | e-value ≤ 1e-5; `-max_target_seqs 5`; 17 curated coral reference sequences |
| pI and MW | BioPython v1.87 ProteinAnalysis | Lehninger pKa scale (default) |
| KEGG pathway names | KEGG REST API | Batch retrieval; 10 KOs per call; 0.2s interval between calls |

### Classification thresholds

| Classification | Threshold / rule |
|---------------|-----------------|
| Gold tier | OrthoFinder ortholog **AND** RBH SwissProt hit (both required) |
| Silver tier | OrthoFinder ortholog **OR** RBH SwissProt hit |
| Bronze tier | Tier 3 HMM hit **OR** informative eggNOG keyword **OR** IPS Pfam hit |
| CellChat Ligand | SignalP class = SP **AND** tmhmm_n_tmrs = 0 |
| CellChat Receptor | In receptor group **OR** (SP + tmhmm_n_tmrs ≥ 1) |
| CellChat Receptor_candidate | tmhmm_n_tmrs ≥ 1, not in known receptor group |

See [`docs/classification_system.md`](docs/classification_system.md) for full classification logic and protein group definitions.


---

## Abbreviations

| Abbreviation | Definition |
|--------------|-----------|
| AMP | Antimicrobial peptide |
| BLAST | Basic Local Alignment Search Tool |
| BSUB | Batch submission command for the LSF scheduler |
| COG | Clusters of Orthologous Groups |
| DBD | DNA-binding domain |
| ECM | Extracellular matrix |
| ESM | Evolutionary Scale Modeling (protein language model used by DeepTMHMM) |
| GFF3 | General Feature Format version 3 |
| GMP | Granulocyte-mast cell precursor |
| GO | Gene Ontology |
| GPCR | G protein-coupled receptor |
| HMM | Hidden Markov Model |
| HPC | High-performance computing cluster |
| HSC | Haematopoietic stem cell |
| HSP | Heat shock protein |
| IPS | InterProScan |
| JID | Job identifier assigned by the LSF scheduler |
| KEGG | Kyoto Encyclopedia of Genes and Genomes |
| KO | KEGG Orthology identifier |
| LSF | Load Sharing Facility (HPC job scheduler) |
| MACPF | Membrane attack complex/perforin domain |
| MSA | Multiple sequence alignment |
| MW | Molecular weight (reported in kDa) |
| NLR | NOD-like receptor |
| OAS | Oligoadenylate synthase |
| OG | Orthogroup (as defined by OrthoFinder) |
| OOM | Out of memory |
| pI | Isoelectric point |
| PRR | Pattern recognition receptor |
| RBH | Reciprocal Best Hit |
| SCRiP | Short Cysteine-Rich Protein (coral-specific venom/defence protein family) |
| SOMP | Skeletal Organic Matrix Protein (coral calcification-associated) |
| SP | Signal peptide (as classified by SignalP 6.0) |
| SRCR | Scavenger receptor cysteine-rich domain |
| STING | Stimulator of Interferon Genes |
| TF | Transcription factor |
| TGF | Transforming growth factor |
| TLR | Toll-like receptor |
| TM | Transmembrane |
| TRAF | TNF receptor-associated factor |
| TRIM | Tripartite motif protein |
| TSV | Tab-separated values |


---

## References

- Stephens TJ, et al. Widespread dark gene evolution in stony corals. *Genome Biol Evol.* 2026;18(4):evag072. https://doi.org/10.1093/gbe/evag072
- Levy S, et al. A single-cell atlas of the coral *Stylophora pistillata*. *Cell.* 2021;184(10):2454–2468. https://doi.org/10.1016/j.cell.2021.04.005
- Helgoe J, et al. Cnidarian innate immunity. *Biol Rev.* 2024. https://doi.org/10.1111/brv.13077
- Lian et al. Coral symbiosis gene regulation. *Sci Adv.* 2025.
- Quigley et al. Coral calcification mechanisms. *Sci Adv.* 2025.

---

## Citation

If you use this pipeline or the annotation data, please cite:

> Wong K. (2026). coral-func-annotation-pipeline. GitHub. https://github.com/kevinhwong1/coral-func-annotation-pipeline

and the framework paper:

> Stephens TJ, et al. Widespread dark gene evolution in stony corals. *Genome Biol Evol.* 2026;18(4):evag072. https://doi.org/10.1093/gbe/evag072

and the individual tool publications listed in the Tool Versions table above.

---

## License

Pipeline scripts: GPL-3.0
Annotation data: CC-BY 4.0
See `LICENSE` for details.
