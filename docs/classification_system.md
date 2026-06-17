# Classification System

This document describes the four classification axes applied to each gene in the master annotation table: confidence tier, annotation status, protein group, and CellChat role.

---

## Contents

- [Confidence Tier](#1-confidence-tier-confidence_tier)
- [Annotation Status](#2-annotation-status-annotation_status)
- [Protein Group](#3-protein-group-protein_group-protein_subgroup-protein_group_secondary)
- [CellChat Role](#4-cellchat-role-cellchat_role)
- [Evidence Priority for best_human_gene_name](#evidence-priority-for-best_human_gene_name)
- [Tool Citations](#tool-citations)
- [Scientific Framework References](#scientific-framework-references)

---

## 1. Confidence Tier (`confidence_tier`)

Reflects the strength of functional evidence from orthology and homology sources. Tiers are mutually exclusive and assigned in priority order.

| Tier | Criteria | Evidence sources |
|------|----------|-----------------|
| **Gold** | OrthoFinder ortholog **AND** RBH SwissProt hit (both required) | `of_human_ortholog` + `rbh_gene_name` |
| **Silver** | OrthoFinder ortholog **OR** RBH SwissProt hit (one of the two) | `of_human_ortholog` or `rbh_gene_name` |
| **Bronze** | Tier 3 HMM domain hit **OR** informative eggNOG keyword **OR** IPS Pfam hit | `tier3_domain_hits`, `eggnog_description`, `interpro_pfam` |
| **Informative** | Tier 4 curated coral family hit only (no ortholog or domain evidence) | `tier4_blast_hit` |
| **Unclassified** | No evidence from any annotation source | |

---

## 2. Annotation Status (`annotation_status`)

Reflects how well a gene is described in public databases, following the framework of Stephens et al. 2026.

| Status | Criteria |
|--------|----------|
| **Light** | Has a named human ortholog OR an informative eggNOG functional description |
| **Partial** | Has a named Pfam domain but no ortholog and no informative description |
| **Dark** | No ortholog, no informative description, and no named domain; gene is functionally uncharacterised |

> Stephens TJ, et al. Widespread dark gene evolution in stony corals. *Genome Biol Evol.* 2026;18(4):evag072. https://doi.org/10.1093/gbe/evag072

---

## 3. Protein Group (`protein_group`, `protein_subgroup`, `protein_group_secondary`)

A functional classification into 24 biologically meaningful groups, designed for use in Seurat metadata and downstream cell communication analyses. Groups are assigned based on a hierarchical evidence cascade: OrthoFinder orthologs → RBH gene names → Tier 3 HMM profile names → eggNOG descriptions → IPS Pfam hits.

### 24 Protein Groups

| Group | Category | Basis for assignment |
|-------|----------|---------------------|
| `GPCR` | Signalling | Pfam: 7tm_1/7tm_2/7tm_3; OG keyword |
| `ECM_Adhesion` | Structural | Pfam: Collagen, FN3, Cadherin, Integrin |
| `Transcription_Factor` | Nuclear | AnimalTFDB DBD HMM hits; TF-related OGs |
| `Kinase` | Signalling | Pfam: Pkinase, HMMR; kinase OGs |
| `Signalling` | Signalling | Wnt, Notch, Hedgehog, TGF-β pathway OGs |
| `Innate_Immunity_PRR` | Immunity | NLR, TLR, SRCR, Pentaxin, lectin-like PRRs |
| `Innate_Immunity_Signalling` | Immunity | MyD88, TRAF, NF-κB, STING pathway |
| `Innate_Immunity_Antiviral` | Immunity | TRIM, DEAD-box helicases, OAS-like |
| `Innate_Immunity_Effector` | Immunity | Complement, perforin, MACPF, ResIII |
| `Lectin_Receptor` | Immunity/Signalling | C-type lectin, galectin domain OGs |
| `Ion_Channel` | Membrane transport | Pfam: Ion_trans, ligand-gated channel OGs |
| `Calcification_Universal` | Calcification | Carbonic anhydrase, bicarbonate transport |
| `Calcification_IonTransport` | Calcification | Ca2+/H+ exchangers, V-ATPase |
| `Calcification_Signalling` | Calcification | Calmodulin, Ca2+-sensing receptors |
| `Symbiosis` | Symbiosis | Autophagy (ATG), lysosomal, sphingolipid |
| `HSP_Chaperone` | Stress response | HSP70, HSP90, HSP40, GroEL/GroES |
| `Neuron_Neuropeptide` | Neuronal | Neuropeptide precursors, synaptic OGs |
| `AMP_Venom` | Defence | Antimicrobial peptide and venom HMMs |
| `Stem_Cell_GMP` | Stem cell | Granulocyte/mast cell precursor markers |
| `Stem_Cell_Proliferation` | Stem cell | Proliferation and cell cycle OGs |
| `Stem_Cell_HSC` | Stem cell | Haematopoietic stem cell-like OGs |
| `Tier4_SCRiP` | Coral-specific | SCRiP (Short Cysteine-Rich Protein) Tier 4 BLAST |
| `Tier4_Neuropep` | Coral-specific | Coral neuropeptide Tier 4 BLAST |
| `Tier4_SOMP` | Coral-specific | Skeletal Organic Matrix Protein Tier 4 BLAST |
| `Unclassified` | — | No group-level evidence |

### Subgroups (`protein_subgroup`)

A finer classification within each group (e.g. `Pentraxin_CRP`, `GPCR_ClassA`, `HSP70`). Populated where sufficient evidence exists; NaN otherwise.

### Dual-Role Genes (`protein_group_secondary`)

A small subset of genes show strong evidence for two functional roles (e.g. a kinase with PRR domain architecture, or a GPCR with neuropeptide precursor features). These receive a primary `protein_group` and a secondary `protein_group_secondary`. The secondary field is NaN for all single-role genes.

---

## 4. CellChat Role (`cellchat_role`)

Classifies each protein's predicted role in intercellular communication, for direct use in CellChat or NicheNet analyses.

| Role | Assignment criteria |
|------|---------------------|
| **Ligand** | SignalP predicts a signal peptide (`SP`) AND `tmhmm_n_tmrs == 0` (secreted, no TM anchor) |
| **Receptor** | Gene is in a receptor protein group (GPCR, Lectin_Receptor, Ion_Channel) OR (SignalP `SP` + `tmhmm_n_tmrs ≥ 1`) |
| **Receptor_candidate** | `tmhmm_n_tmrs ≥ 1`; not in a known receptor group. Predicted membrane-spanning, putative receptor. |
| **Neither** | Cytoplasmic or nuclear proteins (no SP, no TM helices, not in receptor group) |

---

## Evidence Priority for `best_human_gene_name`

When multiple sources provide a human gene name, the following priority order is used:

1. OrthoFinder human ortholog (`of_human_ortholog`): highest confidence
2. RBH SwissProt gene name (`rbh_gene_name`): reciprocal best hit
3. eggNOG human gene (`eggnog_human_gene`): OG-level assignment

---

## Tool Citations

| Tool | Version | Citation |
|------|---------|----------|
| eggNOG-mapper | v2 | Cantalapiedra et al. *Mol Biol Evol.* 2021. https://doi.org/10.1093/molbev/msab293 |
| OrthoFinder | v2 | Emms & Kelly. *Genome Biol.* 2019;20:238. https://doi.org/10.1186/s13059-019-1832-y |
| InterProScan | 5.78-109.0 | Jones et al. *Bioinformatics.* 2014;30:1236–1240. https://doi.org/10.1093/bioinformatics/btu031 |
| DIAMOND | v2.2.1 | Buchfink et al. *Nat Methods.* 2021;18:366–368. https://doi.org/10.1038/s41592-021-01101-x |
| SignalP | 6.0 | Teufel et al. *Nat Biotechnol.* 2022;40:1023–1025. https://doi.org/10.1038/s41587-021-01156-3 |
| DeepTMHMM | 1.0 | Hallgren et al. *bioRxiv.* 2022. https://doi.org/10.1101/2022.04.08.487609 |
| HMMER | 3.4 | Eddy SR. *PLoS Comput Biol.* 2011;7(10):e1002195. https://doi.org/10.1371/journal.pcbi.1002195 |
| BioPython | v1.87 | Cock et al. *Bioinformatics.* 2009;25:1422–1423. https://doi.org/10.1093/bioinformatics/btp163 |

---

## Scientific Framework References

- Stephens TJ, et al. Widespread dark gene evolution in stony corals. *Genome Biol Evol.* 2026;18(4):evag072. https://doi.org/10.1093/gbe/evag072
- Levy S, et al. A single-cell atlas of the coral *Stylophora pistillata*. *Cell.* 2021;184(10):2454–2468. https://doi.org/10.1016/j.cell.2021.04.005
- Helgoe J, et al. Cnidarian immunity. *Biol Rev.* 2024. https://doi.org/10.1111/brv.13077
- Lian et al. Coral symbiosis gene regulation. *Sci Adv.* 2025.
- Quigley et al. Coral calcification mechanisms. *Sci Adv.* 2025.
