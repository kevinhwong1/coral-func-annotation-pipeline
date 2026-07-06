#!/usr/bin/env bash
# ============================================================
# run_config.sh  v4  —  MASTER CONFIGURATION for annotation pipeline
#
# IMPORTANT: This file sets variables ONLY. It does NOT activate
# conda environments — that is the responsibility of each script
# that sources this file.
#
# Source at the top of every pipeline script:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/run_config.sh" "${SPECIES}" "${PROTEOME_PATH:-}"
#
# Usage:
#   source run_config.sh Gfas
#   source run_config.sh Acer /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa
#
# New species that are not listed in the case block can be run by
# passing paths directly via submit_pipeline.sh flags — no edits
# to this file are required.
# ============================================================

# ------------------------------------------------------------
# TOP-LEVEL PATHS
# ------------------------------------------------------------
export NETHOME="/nethome/kxw755"
export BASE="/scratch/dark_genes/annotation_pipeline"
export DB_DIR="${BASE}/databases"

# Conda — NOTE: activation is done explicitly in each script
export CONDA_BASE="${NETHOME}/anaconda3"
export CONDA_ENV_ANNOT="annotation_env"
export CONDA_INIT="source ${NETHOME}/anaconda3/etc/profile.d/conda.sh"

# LSF
export EMAIL="kxw755@earth.miami.edu"
export LSF_PROJECT="dark_genes"

# ------------------------------------------------------------
# DATABASE PATHS
# ------------------------------------------------------------
export SWISSPROT_FASTA="${DB_DIR}/swissprot/uniprot_sprot.fasta"
export HUMAN_FASTA="${DB_DIR}/reference_proteomes/human_swissprot.faa"
export MOUSE_FASTA="${DB_DIR}/reference_proteomes/mouse_swissprot.faa"
export EGGNOG_DB="${DB_DIR}/eggnog"
export EGGNOG_DATA_DIR="${DB_DIR}/eggnog"
export INTERPROSCAN_DIR="${DB_DIR}/interproscan/interproscan-5.78-109.0"
export INTERPROSCAN_BIN="${INTERPROSCAN_DIR}/interproscan.sh"
export ANIMALTFDB_HMM="${DB_DIR}/animaltfdb/all_DBD.hmm"
export TIER3_HMM="${DB_DIR}/tier3_hmms/coral_functional.hmm"
export TIER4_FASTA="${DB_DIR}/tier4_references/tier4_coral_families.fasta"
export NVEC_PROTEOME="${NETHOME}/genomes/NV2_Nvec/NV2g.20240221.protein.fa"

# ------------------------------------------------------------
# TOOL PATHS
# ------------------------------------------------------------
export DEEPTMHMM_DIR="${NETHOME}/DeepTMHMM-Academic-License-v1.0"
export DEEPTMHMM_PREDICT="${DEEPTMHMM_DIR}/predict.py"

# ------------------------------------------------------------
# SPECIES CONFIGURATION
# Call: configure_species SPECIES [PROTEOME_PATH]
# Sets all run-specific path variables.
# Does NOT activate conda or run any commands.
#
# To add a new species, add a case block below with:
#   QUERY_PROTEOME  — path to protein FASTA
#   SPECIES_GFF3    — path to GFF3 annotation file
#   GENOME_VERSION  — genome version string for annotation summary
#   PROTEOME_SOURCE — source URL for annotation summary
#
# All four are optional if passed via submit_pipeline.sh flags.
# ------------------------------------------------------------
configure_species() {
    local SPECIES_ARG="${1:-Gfas}"
    local PROTEOME_ARG="${2:-}"

    export SPECIES="${SPECIES_ARG}"

    # Safe defaults for new variables — prevents set -u failures
    export SPECIES_GFF3="${SPECIES_GFF3:-}"
    export GENOME_VERSION="${GENOME_VERSION:-}"
    export PROTEOME_SOURCE="${PROTEOME_SOURCE:-}"
    export STAGING_DIR="${STAGING_DIR:-${NETHOME}/github_upload}"
    export GENE_TO_PROTEIN="${GENE_TO_PROTEIN:-}"

    # Default proteome paths per species — add new species here
    case "${SPECIES}" in
        Gfas)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/Gfas_v1/gfas_1.0.proteins.fasta}"
            export SPECIES_GFF3="${SPECIES_GFF3:-/nethome/kxw755/genomes/Gfas_v1/gfas_1.0.genes.gff3}"
            export GENOME_VERSION="${GENOME_VERSION:-Gfas_v1.0}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-http://gfas.reefgenomics.org/}"
            ;;
        Pdam)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/pdam_genome/pdam_proteins.fasta}"
            export SPECIES_GFF3="${SPECIES_GFF3:-/nethome/kxw755/genomes/pdam_genome/pdam_annotation.gff3}"
            export GENOME_VERSION="${GENOME_VERSION:-Pdam_v1.0}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-http://pdam.reefgenomics.org/}"
            ;;
        Acer)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/GCA_032359415.1_Acer/acer_protein.faa}"
            export SPECIES_GFF3="${SPECIES_GFF3:-${NETHOME}/genomes/GCA_032359415.1_Acer/genomic.gff}"
            export GENOME_VERSION="${GENOME_VERSION:-GCA_032359415.1}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/}"
            export GENE_TO_PROTEIN="${GENE_TO_PROTEIN:-${NETHOME}/genomes/GCA_032359415.1_Acer/acer_gene_to_protein_combined.csv}"
            ;;
        Ofav)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/longreads_um_ofav_complete/Orbicella_faveolata_gen_17.proteins.fa}"
            export SPECIES_GFF3="${SPECIES_GFF3:-${NETHOME}/genomes/longreads_um_ofav_complete/Orbicella_faveolata_gen_17.gff3}"
            export GENOME_VERSION="${GENOME_VERSION:-Ofav_gen_17}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-University of Miami long-read assembly (internal)}"
            ;;
        Amur)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/GCF_036669905.1_Amur/amur_protein.faa}"
            export SPECIES_GFF3="${SPECIES_GFF3:-${NETHOME}/genomes/GCF_036669905.1_Amur/genomic.gff}"
            export GENOME_VERSION="${GENOME_VERSION:-GCF_036669905.1}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669905.1/}"
            export GENE_TO_PROTEIN="${GENE_TO_PROTEIN:-${NETHOME}/genomes/GCF_036669905.1_Amur/amur_gene_to_protein.csv}"
            ;;
        Amil)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/Amil/Amil_long.pep.fasta}"
            export GENOME_VERSION="${GENOME_VERSION:-Amil_v2.1}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://github.com/sebepedroslab/oculina-coral-sc-atlas}"
            ;;
        Spis)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/Spis/Spis_long.pep.fasta}"
            export GENOME_VERSION="${GENOME_VERSION:-Spis_v1.0}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://github.com/sebepedroslab/oculina-coral-sc-atlas}"
            ;;
        Opat)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/Opat/Ocupat_long.pep.fasta}"
            export GENOME_VERSION="${GENOME_VERSION:-Opat_v1.0}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://github.com/sebepedroslab/oculina-coral-sc-atlas}"
            ;;
        Xspp)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NETHOME}/genomes/Xspp/xenSp1.proteins.fa}"
            export GENOME_VERSION="${GENOME_VERSION:-xenSp1}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://cmo.carnegiescience.edu/endosymbiosis/genome/}"
            ;;
        Nvec)
            export QUERY_PROTEOME="${PROTEOME_ARG:-${NVEC_PROTEOME}}"
            export SPECIES_GFF3="${SPECIES_GFF3:-${NETHOME}/genomes/NV2_Nvec/NV2g.20240221.gff}"
            export GENOME_VERSION="${GENOME_VERSION:-NV2g.20240221}"
            export PROTEOME_SOURCE="${PROTEOME_SOURCE:-https://simrbase.stowers.org/nematostella}"
            ;;
        *)
            if [[ -z "${PROTEOME_ARG}" ]]; then
                echo "ERROR: Unknown species '${SPECIES}'. Supply proteome path via --proteome flag." >&2
                echo "Usage: bash submit_pipeline.sh ${SPECIES} --proteome /path/to/proteome.fasta" >&2
                return 1
            fi
            export QUERY_PROTEOME="${PROTEOME_ARG}"
            ;;
    esac

    # Validate proteome exists
    if [[ ! -f "${QUERY_PROTEOME}" ]]; then
        echo "ERROR: Proteome not found: ${QUERY_PROTEOME}" >&2
        return 1
    fi

    # Dated run directory
    export RUNDATE=$(date +%Y%m%d)
    export RUN_DIR="${BASE}/runs/${SPECIES}_${RUNDATE}"
    export LOG_DIR="${RUN_DIR}/logs"
    export SCRIPT_DIR_ANNOT="${BASE}/scripts"

    # Per-step output directories
    export EGGNOG_OUT="${RUN_DIR}/01_eggnog"
    export OF_OUT="${RUN_DIR}/02_orthofinder"
    export IPS_OUT="${RUN_DIR}/03_interproscan"
    export RBH_OUT="${RUN_DIR}/04_rbh_swissprot"
    export SIGNALP_OUT="${RUN_DIR}/05_signalp"
    export TMHMM_OUT="${RUN_DIR}/05_tmhmm"
    export TIER3_OUT="${RUN_DIR}/06_tier3_hmmsearch"
    export TIER4_OUT="${RUN_DIR}/07_tier4_blast"
    export FINAL_OUT="${RUN_DIR}/08_final"

    # OrthoFinder input directory
    export OF_INPUT_DIR="${RUN_DIR}/02_orthofinder/input_proteomes"

    # Key output files used by merge step
    export EGGNOG_ANNOT="${EGGNOG_OUT}/${SPECIES}.emapper.annotations"
    export OF_RESULTS_DIR="${OF_OUT}/Results_OrthoFinder"
    export IPS_TSV="${IPS_OUT}/${SPECIES}_interproscan.tsv"
    export RBH_TSV="${RBH_OUT}/${SPECIES}_rbh_hits.tsv"
    export SIGNALP_TSV="${SIGNALP_OUT}/prediction_results.txt"
    export TMHMM_GFF3="${TMHMM_OUT}/TMRs.gff3"
    export TIER3_DOMTBL="${TIER3_OUT}/${SPECIES}_tier3_domains.domtblout"
    export TIER4_BLAST="${TIER4_OUT}/${SPECIES}_tier4_blast.tsv"
    export MASTER_TABLE="${FINAL_OUT}/${SPECIES}_master_annotation.tsv"
    export GENE_LISTS_DIR="${FINAL_OUT}/gene_lists"
}

# Auto-configure if called with arguments
if [[ -n "${1:-}" ]]; then
    configure_species "${1}" "${2:-}"
fi
