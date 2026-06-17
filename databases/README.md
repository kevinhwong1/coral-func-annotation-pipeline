# Databases

## Contents

- [Database List](#database-list)
- [Notes](#notes)

---

## Database List

Database files are not hosted in this repository due to size. Download from the sources below before running the pipeline. Custom databases are included in the repo and linked below.

| Database | Version | Size | Download source | Pipeline path |
|----------|---------|------|----------------|---------------|
| eggNOG 5.0 | 5.0 | ~39 GB | http://eggnog-mapper.embl.de | `databases/eggnog/` |
| UniProt Swiss-Prot | Jan 2026 | ~275 MB | https://www.uniprot.org/downloads | `databases/swissprot/uniprot_sprot.fasta` |
| InterProScan | 5.78-109.0 | ~20 GB | https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/ | `databases/interproscan/` |
| AnimalTFDB DBD HMMs | 4.0 | <1 MB | https://guolab.wchscu.cn/AnimalTFDB4/ | `databases/animaltfdb/all_DBD.hmm` |
| Human SwissProt proteome | Jan 2026 | ~3 MB | https://www.uniprot.org/proteomes/UP000005640 | `databases/reference_proteomes/human_swissprot.faa` |
| Mouse SwissProt proteome | Jan 2026 | ~3 MB | https://www.uniprot.org/proteomes/UP000000589 | `databases/reference_proteomes/mouse_swissprot.faa` |
| Tier 3 HMMs (custom) | — | <1 MB | [databases/tier3_hmms/coral_functional.hmm](tier3_hmms/coral_functional.hmm) | `databases/tier3_hmms/coral_functional.hmm` |
| Tier 4 reference sequences (custom) | — | <1 MB | [databases/tier4_references/tier4_coral_families.fasta](tier4_references/tier4_coral_families.fasta) | `databases/tier4_references/tier4_coral_families.fasta` |
| Pfam name mappings (custom) | — | <1 MB | [databases/pfam_names.tsv](pfam_names.tsv) | `databases/pfam_names.tsv` |

## Notes

- **Tier 3 HMMs**: 179 coral/metazoan functional profiles covering immunity, calcification, GPCRs, kinases, ECM, stem cell, and stress response domains. Compiled from Pfam and AnimalTFDB4. See [`docs/classification_system.md`](../docs/classification_system.md) for the full profile list and group assignments.
- **Tier 4 references**: 17 curated coral-specific sequences (SCRiPs, SOMPs, neuropeptides, cnidocyte proteins) downloaded from UniProt. See [`docs/classification_system.md`](../docs/classification_system.md) for accessions and sources.
- **Pfam name mappings**: 6,499 Pfam accession → human-readable name mappings, extracted from the first InterProScan run output. Regenerate with the one-liner in `pipeline/00_setup_annotation_env.sh` after running InterProScan on a new species.
- **SwissProt**: downloaded from `ftp.uniprot.org/pub/databases/uniprot/current_release/` (no version pinned; date of download recorded above).
