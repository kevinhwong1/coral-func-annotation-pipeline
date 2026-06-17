#!/usr/bin/env bash
# ============================================================
# 00_setup.sh
#
# ONE-TIME setup on Pegasus login node.
# Run interactively (NOT as a BSUB job):
#   bash scripts/00_setup.sh
#
# What this does:
#   1. Creates the full directory skeleton on scratch
#   2. Creates the annotation conda environment
#   3. Downloads all databases (takes 3-6 hrs, network dependent)
#   4. Compiles the Tier 3 HMM database from individual Pfam profiles
#   5. Builds the Tier 4 curated coral sequence FASTA
#   6. Validates everything is in place
#
# IMPORTANT: Only needs to run ONCE. All future species runs
#            reuse the same databases.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/run_config.sh"

echo "============================================================"
echo " Annotation Pipeline: One-Time Setup"
echo " Base directory: ${BASE}"
echo "============================================================"
date

# ============================================================
# 1. Create directory skeleton
# ============================================================
echo ""
echo "--- Creating directory skeleton ---"

mkdir -p "${BASE}/databases/swissprot"
mkdir -p "${BASE}/databases/reference_proteomes"
mkdir -p "${BASE}/databases/eggnog"
mkdir -p "${BASE}/databases/interproscan"
mkdir -p "${BASE}/databases/animaltfdb"
mkdir -p "${BASE}/databases/tier3_hmms"
mkdir -p "${BASE}/databases/tier4_references"
mkdir -p "${BASE}/runs"
mkdir -p "${BASE}/logs"

echo "  Directory skeleton created."

# ============================================================
# 2. Create conda environment
# ============================================================
echo ""
echo "--- Setting up conda environment: ${CONDA_ENV_ANNOT} ---"

source ~/anaconda3/bin/activate

if conda env list | grep -q "^${CONDA_ENV_ANNOT} "; then
    echo "  Environment '${CONDA_ENV_ANNOT}' already exists — skipping creation."
    echo "  To rebuild: conda env remove -n ${CONDA_ENV_ANNOT}"
else
    echo "  Creating conda environment: ${CONDA_ENV_ANNOT}"
    conda create -n "${CONDA_ENV_ANNOT}" python=3.11 -y

    conda activate "${CONDA_ENV_ANNOT}"

    # Bioinformatics tools
    conda install -n "${CONDA_ENV_ANNOT}" -c bioconda -c conda-forge \
        hmmer \
        blast \
        diamond \
        -y

    # Python packages for merge/annotate step
    pip install \
        pandas \
        numpy \
        biopython \
        openpyxl \
        tqdm \
        --quiet

    echo "  Conda environment created successfully."
fi

conda activate "${CONDA_ENV_ANNOT}"

# ============================================================
# 3. Download UniProt SwissProt (reviewed proteins only)
# ============================================================
echo ""
echo "--- Downloading UniProt SwissProt ---"
echo "  Target: ${SWISSPROT_FASTA}"

if [[ -s "${SWISSPROT_FASTA}" ]]; then
    echo "  Already exists ($(du -sh ${SWISSPROT_FASTA} | cut -f1)). Skipping."
else
    echo "  Downloading... (~250 MB)"
    wget -q --show-progress \
        -O "${SWISSPROT_FASTA}.gz" \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip "${SWISSPROT_FASTA}.gz"
    echo "  Done: $(grep -c '^>' ${SWISSPROT_FASTA}) sequences"
fi

# ============================================================
# 4. Extract human + mouse SwissProt subsets for OrthoFinder
# ============================================================
echo ""
echo "--- Extracting human and mouse SwissProt subsets ---"

python3 - <<'PY'
import os, re, sys
from pathlib import Path

swissprot = os.environ["SWISSPROT_FASTA"]
human_out  = os.environ["HUMAN_FASTA"]
mouse_out  = os.environ["MOUSE_FASTA"]

def extract_organism(sprot_fasta, out_fasta, os_code):
    """Extract sequences for one organism code (e.g. HUMAN, MOUSE)."""
    # SwissProt headers: >sp|P08107|HSPA1_HUMAN Heat shock 70 kDa...
    # We match the entry name suffix (e.g. _HUMAN)
    pattern = re.compile(rf'\|[A-Z0-9]+_{os_code}\b')
    n = 0
    with open(sprot_fasta) as f_in, open(out_fasta, 'w') as f_out:
        write = False
        for line in f_in:
            if line.startswith('>'):
                write = bool(pattern.search(line))
                if write:
                    # Rename header to just gene name for OrthoFinder clarity
                    # e.g. >sp|P08107|HSPA1_HUMAN ... → >HSPA1_HUMAN P08107
                    m = re.match(r'>sp\|([A-Z0-9]+)\|(\S+)', line)
                    if m:
                        acc, entry = m.group(1), m.group(2)
                        line = f'>{entry} {acc}\n'
                        n += 1
            if write:
                f_out.write(line)
    print(f"  {os_code}: {n:,} sequences → {out_fasta}")

for out_path in [human_out, mouse_out]:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

if not Path(human_out).exists() or os.path.getsize(human_out) == 0:
    extract_organism(swissprot, human_out, "HUMAN")
else:
    print(f"  HUMAN: already exists, skipping")

if not Path(mouse_out).exists() or os.path.getsize(mouse_out) == 0:
    extract_organism(swissprot, mouse_out, "MOUSE")
else:
    print(f"  MOUSE: already exists, skipping")
PY

# ============================================================
# 5. Make BLAST database for SwissProt (for RBH step)
# ============================================================
echo ""
echo "--- Making SwissProt BLAST database ---"

if [[ -f "${SWISSPROT_FASTA}.pdb" ]]; then
    echo "  BLAST db already exists. Skipping."
else
    makeblastdb \
        -in "${SWISSPROT_FASTA}" \
        -dbtype prot \
        -out "${SWISSPROT_FASTA}" \
        -parse_seqids \
        -title "UniProt_SwissProt_reviewed"
    echo "  BLAST database created."
fi

# ============================================================
# 6. Make BLAST database for query proteome (for RBH reverse)
# Note: This will need to be run again for each new species.
# The submit_pipeline.sh script handles this per-run.
# ============================================================
echo ""
echo "--- Making BLAST database for Gfas proteome ---"
makeblastdb \
    -in "${QUERY_PROTEOME}" \
    -dbtype prot \
    -out "${QUERY_PROTEOME}" \
    -parse_seqids \
    -title "${SPECIES}_proteome" \
    2>/dev/null || echo "  NOTE: BLAST db for query proteome will be made per-run."

# ============================================================
# 7. Download eggNOG database
# ============================================================
echo ""
echo "--- Downloading eggNOG database ---"
echo "  Target: ${EGGNOG_DB}"
echo "  WARNING: This is ~50 GB and will take 1-3 hours."

mkdir -p "${EGGNOG_DB}"

if [[ -f "${EGGNOG_DB}/eggnog.db" ]]; then
    echo "  eggnog.db already exists ($(du -sh ${EGGNOG_DB}/eggnog.db | cut -f1)). Skipping."
else
    conda activate "${CONDA_ENV_ANNOT}"
    download_eggnog_data.py \
        --data_dir "${EGGNOG_DB}" \
        -y \
        2>&1 | tee "${BASE}/logs/eggnog_download.log"
    echo "  eggNOG database download complete."
fi

# ============================================================
# 8. Download and install InterProScan
# ============================================================
echo ""
echo "--- Setting up InterProScan ---"
echo "  Target: ${INTERPROSCAN_DIR}"
echo "  WARNING: This is ~20 GB and will take 1-2 hours."

IPS_VERSION="5.71-103.0"
IPS_TARBALL="${DB_DIR}/interproscan/interproscan-${IPS_VERSION}-64-bit.tar.gz"

if [[ -f "${INTERPROSCAN_BIN}" ]]; then
    echo "  InterProScan already installed. Skipping."
else
    mkdir -p "${DB_DIR}/interproscan"
    cd "${DB_DIR}/interproscan"

    echo "  Downloading InterProScan ${IPS_VERSION}..."
    wget -q --show-progress \
        "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPS_VERSION}/interproscan-${IPS_VERSION}-64-bit.tar.gz" \
        -O "${IPS_TARBALL}"

    echo "  Extracting..."
    tar -xzf "${IPS_TARBALL}"
    rm "${IPS_TARBALL}"

    # Java is required — check for system Java first
    if ! java -version 2>&1 | grep -q "version"; then
        echo "  Java not found in PATH. Installing via conda..."
        conda install -n "${CONDA_ENV_ANNOT}" -c conda-forge openjdk=17 -y
    fi

    echo "  InterProScan installed: ${INTERPROSCAN_BIN}"
fi

# ============================================================
# 9. Download AnimalTFDB4 HMM profiles
# ============================================================
echo ""
echo "--- Downloading AnimalTFDB4 HMM profiles ---"

TFDB_HMM_DIR="${DB_DIR}/animaltfdb/hmm_profiles"
mkdir -p "${TFDB_HMM_DIR}"

if [[ -f "${ANIMALTFDB_HMM}" ]]; then
    echo "  AnimalTFDB HMMs already compiled. Skipping."
else
    echo "  Downloading from AnimalTFDB4..."
    # Primary mirror
    wget -q --show-progress \
        -O "${TFDB_HMM_DIR}/HMM_DBD.tar.gz" \
        "https://guolab.wchscu.cn/AnimalTFDB4/static/download/HMM_DBD.tar.gz" \
    || wget -q --show-progress \
        -O "${TFDB_HMM_DIR}/HMM_DBD.tar.gz" \
        "http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/HMM_DBD.tar.gz"

    tar -xzf "${TFDB_HMM_DIR}/HMM_DBD.tar.gz" -C "${TFDB_HMM_DIR}/"

    echo "  Concatenating all TF family HMMs..."
    cat "${TFDB_HMM_DIR}"/*.hmm > "${ANIMALTFDB_HMM}"

    echo "  Pressing (indexing) combined HMM database..."
    hmmpress "${ANIMALTFDB_HMM}"

    echo "  AnimalTFDB HMMs ready: $(grep -c '^NAME' ${ANIMALTFDB_HMM}) profiles"
fi

# ============================================================
# 10. Download Tier 3 Pfam HMMs and compile coral_functional.hmm
# ============================================================
echo ""
echo "--- Compiling Tier 3 coral functional HMM database ---"

PFAM_DIR="${DB_DIR}/tier3_hmms/pfam_individual"
mkdir -p "${PFAM_DIR}"

if [[ -f "${TIER3_HMM}" ]]; then
    echo "  Tier 3 HMM already compiled. Skipping."
else
    echo "  Downloading individual Pfam HMM profiles..."

    # Function to download one Pfam HMM
    download_pfam_hmm() {
        local pfam_id="$1"
        local label="$2"
        local out_file="${PFAM_DIR}/${pfam_id}_${label}.hmm"
        if [[ ! -f "${out_file}" ]]; then
            wget -q \
                "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${pfam_id}/?annotation=hmm" \
                -O "${out_file}" \
            || echo "  WARNING: Could not download ${pfam_id} (${label})"
        fi
    }

    # ── IMMUNITY ────────────────────────────────────────────────
    echo "  Immunity HMMs..."
    download_pfam_hmm "PF01582" "TIR_domain_TLR"
    download_pfam_hmm "PF00560" "LRR_TLR_extracellular"
    download_pfam_hmm "PF05729" "NACHT_NLR"
    download_pfam_hmm "PF13855" "LRR_8_NLR"
    download_pfam_hmm "PF00530" "SRCR_scavenger_receptor"
    download_pfam_hmm "PF00059" "CTLD_C_type_lectin"
    download_pfam_hmm "PF00337" "Galectin"
    download_pfam_hmm "PF15009" "STING_TM"
    download_pfam_hmm "PF00656" "Caspase_domain"
    download_pfam_hmm "PF01823" "MACPF_perforin"
    download_pfam_hmm "PF00275" "Endonuclease_G"
    download_pfam_hmm "PF00270" "DEAD_helicase_RIG_I_Vasa"
    download_pfam_hmm "PF00071" "Ras_GTPase_Rab"
    download_pfam_hmm "PF02022" "Integrase_zinc_TRIM"
    download_pfam_hmm "PF00264" "Tyrosinase_melanin"
    download_pfam_hmm "PF01476" "LysM_lectin_ficolin"

    # ── STEM CELLS ──────────────────────────────────────────────
    echo "  Stem cell HMMs..."
    download_pfam_hmm "PF02171" "Piwi_Argonaute"
    download_pfam_hmm "PF00270" "DEAD_helicase_Vasa_DDX4"
    download_pfam_hmm "PF03417" "Nanos_zinc_finger"
    download_pfam_hmm "PF00806" "Pumilio_PUF"
    download_pfam_hmm "PF00705" "PCNA_proliferation"
    download_pfam_hmm "PF02042" "RVP_replication"

    # ── CALCIFICATION ───────────────────────────────────────────
    echo "  Calcification HMMs..."
    download_pfam_hmm "PF00194" "Carbonic_anhydrase_alpha"
    download_pfam_hmm "PF07836" "Carbonic_anhydrase_beta"
    download_pfam_hmm "PF07565" "SLC4_bicarbonate_transporter"
    download_pfam_hmm "PF01699" "SLC26_Na_bicarbonate"
    download_pfam_hmm "PF00122" "Ca_ATPase_SERCA"
    download_pfam_hmm "PF00137" "V_ATPase_H_pump"
    download_pfam_hmm "PF00093" "Von_Willebrand_D_domain"
    download_pfam_hmm "PF00187" "Chitin_binding_domain"
    download_pfam_hmm "PF00704" "Chitinase"
    download_pfam_hmm "PF01644" "Chitin_synthase"
    download_pfam_hmm "PF00008" "EGF_like_domain"
    download_pfam_hmm "PF00041" "Fibronectin_type_III"
    download_pfam_hmm "PF01389" "Collagen_fibrillar"

    # ── SIGNALLING / RECEPTORS ───────────────────────────────────
    echo "  Signalling / receptor HMMs..."
    download_pfam_hmm "PF00069" "Kinase_domain"
    download_pfam_hmm "PF07714" "Kinase_Tyr_RTK"
    download_pfam_hmm "PF00001" "GPCR_7tm_1"
    download_pfam_hmm "PF10320" "GPCR_7tm_3"
    download_pfam_hmm "PF00520" "Voltage_gated_ion_channel"
    download_pfam_hmm "PF00076" "RRM_RNA_binding_ELAV"

    # ── HSP / CHAPERONES ────────────────────────────────────────
    echo "  HSP/chaperone HMMs..."
    download_pfam_hmm "PF00012" "HSP70_ATPase"
    download_pfam_hmm "PF02976" "HSP70_substrate_binding"
    download_pfam_hmm "PF00118" "HSP60_GroEL"
    download_pfam_hmm "PF00226" "DnaJ_HSP40"
    download_pfam_hmm "PF02518" "HATPase_HSP90"

    # Concatenate all downloaded HMMs
    echo "  Concatenating Tier 3 HMMs..."
    cat "${PFAM_DIR}"/*.hmm > "${TIER3_HMM}"

    echo "  Pressing Tier 3 HMM database..."
    hmmpress "${TIER3_HMM}"

    N_PROFILES=$(grep -c '^NAME' "${TIER3_HMM}")
    echo "  Tier 3 HMM database ready: ${N_PROFILES} profiles"
fi

# ============================================================
# 11. Build Tier 4 curated coral sequence FASTA
# ============================================================
echo ""
echo "--- Building Tier 4 curated coral sequence FASTA ---"

if [[ -f "${TIER4_FASTA}" ]]; then
    N_T4=$(grep -c '^>' "${TIER4_FASTA}")
    echo "  Already exists: ${N_T4} sequences. Skipping."
else
    echo "  Downloading Tier 4 reference sequences from NCBI..."

    python3 - <<'PY'
import os, subprocess, textwrap
from pathlib import Path

tier4_out = os.environ["TIER4_FASTA"]
Path(tier4_out).parent.mkdir(parents=True, exist_ok=True)

# Each entry: (accession, group, species_code, gene_name)
# Header format: >GROUP_Species_GeneName_Accession
SEQUENCES = [
    # ── GALAXIN / SOMP Calcification ─────────────────────────────
    ("Q8I6S1",   "SOMP_Calc", "Gfas", "Galaxin"),            # Galaxea fascicularis galaxin ORIGINAL
    ("Q5I3W9",   "SOMP_Calc", "Amil", "Galaxin"),            # Acropora millepora galaxin
    ("Q5I3W8",   "SOMP_Calc", "Amil", "GalaxinLike1"),       # A. millepora galaxin-like 1
    ("Q5I3W7",   "SOMP_Calc", "Amil", "GalaxinLike2"),       # A. millepora galaxin-like 2
    # CARPs from Stylophora pistillata (Drake et al. 2013)
    ("P86858",   "SOMP_Calc", "Spis", "CARP1"),
    ("P86859",   "SOMP_Calc", "Spis", "CARP2"),
    ("P86860",   "SOMP_Calc", "Spis", "CARP3"),
    ("P86861",   "SOMP_Calc", "Spis", "CARP4"),
    # Additional Spis SOMPs from Drake et al. 2013 / Ramos-Silva et al. 2013
    ("K7PXU5",   "SOMP_Calc", "Spis", "SOM_P1"),
    ("K7Q056",   "SOMP_Calc", "Spis", "SOM_P5"),
    ("K7Q1U8",   "SOMP_Calc", "Spis", "SOM_P14"),
    # ── SCRiPs ───────────────────────────────────────────────────
    ("B7SPD1",   "SCRiP_Venom",  "Mfav", "SCRiP1"),          # Orbicella faveolata
    ("B7SPD2",   "SCRiP_Venom",  "Mfav", "SCRiP2"),
    ("B7SPD3",   "SCRiP_Venom",  "Mfav", "SCRiP3"),
    ("B7X8E0",   "SCRiP_Venom",  "Amil", "SCRiP1"),          # Acropora millepora
    ("B7X8E1",   "SCRiP_Venom",  "Amil", "SCRiP2"),
    ("B7X8E2",   "SCRiP_Venom",  "Amil", "SCRiP3"),
    # ── RFamide neuropeptides (Hexacorallia) ──────────────────────
    ("A7RZI3",   "Neuropep_NP",  "Nvec", "RFamide_prepro"),  # Nematostella vectensis
    # ── GLWamide neuropeptides ────────────────────────────────────
    ("Q86LQ0",   "Neuropep_NP",  "Hm",   "GLWamide_prepro"), # Hydra magnipapillata
    # ── Symbiosis: LePin (lectin-kazal) ──────────────────────────
    ("A0A8B8P0I1","Symbiosis",   "Xesp", "LePin"),           # Xenia sp. LePin
    # ── Minicollagen (cnidocyte) ──────────────────────────────────
    ("P31835",   "Cnidocyte",    "Hvul", "Minicollagen1"),    # Hydra vulgaris
    ("Q9Y045",   "Cnidocyte",    "Nvec", "Minicollagen1"),    # Nematostella
]

# Download sequences from UniProt
accessions = [s[0] for s in SEQUENCES]
label_map  = {s[0]: f">{s[1]}_{s[2]}_{s[3]}_{s[0]}" for s in SEQUENCES}

print(f"  Downloading {len(accessions)} Tier 4 sequences from UniProt...")

import urllib.request, time

with open(tier4_out, "w") as f_out:
    for acc, group, sp, gene in SEQUENCES:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                fasta_text = resp.read().decode("utf-8")
            # Replace header with our standard format
            lines = fasta_text.strip().split("\n")
            new_header = f">{group}_{sp}_{gene}_{acc}"
            f_out.write(new_header + "\n")
            f_out.write("\n".join(lines[1:]) + "\n")
            print(f"    OK: {acc} → {group}_{sp}_{gene}")
        except Exception as e:
            print(f"    WARNING: Could not download {acc}: {e}")
        time.sleep(0.3)  # be polite to UniProt

# Count
with open(tier4_out) as f:
    n = sum(1 for line in f if line.startswith(">"))
print(f"\n  Tier 4 FASTA: {n} sequences written to {tier4_out}")
PY

    # Make BLAST database for Tier 4
    makeblastdb \
        -in "${TIER4_FASTA}" \
        -dbtype prot \
        -out "${TIER4_FASTA}" \
        -parse_seqids \
        -title "Tier4_Coral_Families"

    echo "  Tier 4 BLAST database created."
fi

# ============================================================
# 12. Validate everything
# ============================================================
echo ""
echo "--- Validating setup ---"

MISSING=0

check_file() {
    if [[ -s "$1" ]]; then
        echo "  ✓ $2"
    else
        echo "  ✗ MISSING: $2 ($1)"
        MISSING=$((MISSING + 1))
    fi
}

check_file "${SWISSPROT_FASTA}"    "SwissProt FASTA"
check_file "${HUMAN_FASTA}"        "Human proteome (SwissProt)"
check_file "${MOUSE_FASTA}"        "Mouse proteome (SwissProt)"
check_file "${EGGNOG_DB}/eggnog.db" "eggNOG database"
check_file "${INTERPROSCAN_BIN}"   "InterProScan binary"
check_file "${ANIMALTFDB_HMM}"     "AnimalTFDB HMM database"
check_file "${TIER3_HMM}"          "Tier 3 coral functional HMM"
check_file "${TIER4_FASTA}"        "Tier 4 coral sequences"
check_file "${NVEC_PROTEOME}"      "Nematostella proteome"

if [[ "${MISSING}" -gt 0 ]]; then
    echo ""
    echo "  WARNING: ${MISSING} items missing. Check logs above."
else
    echo ""
    echo "  All databases present and accounted for."
fi

echo ""
echo "============================================================"
echo " 00_setup.sh COMPLETE"
echo " Databases in: ${DB_DIR}"
echo " Ready to run: bash scripts/submit_pipeline.sh Gfas"
echo "============================================================"
date
