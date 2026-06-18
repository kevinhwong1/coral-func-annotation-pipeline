#!/usr/bin/env bash
# ============================================================
# submit_pipeline.sh  v4
#
# Master submission script. Submits all annotation jobs in
# parallel, then submits the merge job with LSF dependencies.
#
# Jobs:
#   01  eggNOG-mapper
#   02  OrthoFinder (MSA mode)
#   03  InterProScan 5.78
#   04  RBH vs SwissProt
#   05  SignalP 6.0 (secretome)
#   05b DeepTMHMM (TM topology) — auto-chunks proteome
#   06  Tier 3 HMM search (179 coral profiles)
#   07  Tier 4 BLASTp (curated coral sequences)
#   08  Merge + annotate (waits on all above)
#
# Usage:
#   bash scripts/submit_pipeline.sh <SPECIES> [options]
#
# Required:
#   $1  SPECIES code (e.g. Gfas, Acer, Spis)
#
# Options (override run_config.sh defaults):
#   --proteome        Path to proteome FASTA
#   --gff3            Path to GFF3 annotation file
#   --genome_version  Genome version string for annotation summary
#   --proteome_source Proteome source URL for annotation summary
#   --staging_dir     Root staging directory for GitHub upload
#
# Example:
#   bash scripts/submit_pipeline.sh Acer \
#       --proteome /nethome/kxw755/genomes/GCA_032359415.1_Acer/acer_protein.faa \
#       --gff3     /nethome/kxw755/genomes/GCA_032359415.1_Acer/genomic.gff \
#       --genome_version  "GCA_032359415.1" \
#       --proteome_source "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/" \
#       --staging_dir     "/nethome/kxw755/github_upload"
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parse arguments ──────────────────────────────────────────
SPECIES="${1:?ERROR: Supply species code. Usage: bash submit_pipeline.sh Gfas}"
shift

# Optional flag overrides
OPT_PROTEOME=""
OPT_GFF3=""
OPT_GENOME_VERSION=""
OPT_PROTEOME_SOURCE=""
OPT_STAGING_DIR=""
OPT_GENE_TO_PROTEIN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --proteome)        OPT_PROTEOME="$2";        shift 2 ;;
        --gff3)            OPT_GFF3="$2";            shift 2 ;;
        --genome_version)  OPT_GENOME_VERSION="$2";  shift 2 ;;
        --proteome_source) OPT_PROTEOME_SOURCE="$2"; shift 2 ;;
        --staging_dir)     OPT_STAGING_DIR="$2";     shift 2 ;;
        --gene_to_protein) OPT_GENE_TO_PROTEIN="$2"; shift 2 ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Load species config ──────────────────────────────────────
source "${SCRIPT_DIR}/run_config.sh" "${SPECIES}" "${OPT_PROTEOME}"

# Apply flag overrides (take priority over run_config.sh)
[[ -n "${OPT_PROTEOME}" ]]        && export QUERY_PROTEOME="${OPT_PROTEOME}"
[[ -n "${OPT_GFF3}" ]]            && export SPECIES_GFF3="${OPT_GFF3}"
[[ -n "${OPT_GENOME_VERSION}" ]]  && export GENOME_VERSION="${OPT_GENOME_VERSION}"
[[ -n "${OPT_PROTEOME_SOURCE}" ]] && export PROTEOME_SOURCE="${OPT_PROTEOME_SOURCE}"
[[ -n "${OPT_STAGING_DIR}" ]]     && export STAGING_DIR="${OPT_STAGING_DIR}"
[[ -n "${OPT_GENE_TO_PROTEIN}" ]] && export GENE_TO_PROTEIN="${OPT_GENE_TO_PROTEIN}"

echo "============================================================"
echo " Annotation Pipeline v4 Submission"
echo " Species:         ${SPECIES}"
echo " Proteome:        ${QUERY_PROTEOME}"
echo " GFF3:            ${SPECIES_GFF3:-not provided}"
echo " Genome version:  ${GENOME_VERSION:-not provided}"
echo " Staging dir:     ${STAGING_DIR:-not set}"
echo " Run dir:         ${RUN_DIR}"
echo "============================================================"
date

# ── Create run directory structure ──────────────────────────
mkdir -p "${RUN_DIR}"/{logs,01_eggnog,02_orthofinder,03_interproscan,\
04_rbh_swissprot,05_signalp,05_tmhmm/chunks,05_tmhmm/results,\
06_tier3_hmmsearch,07_tier4_blast,08_final/gene_lists,py_helpers}

# Save scripts into run dir for reproducibility
cp "${SCRIPT_DIR}/run_config.sh"        "${RUN_DIR}/"

# Record latest run
echo "${SPECIES}_${RUNDATE}" > "${BASE}/runs/LATEST_RUN.txt"

# ── Activate conda ───────────────────────────────────────────
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env

# ── Make BLAST database for query proteome ──────────────────
echo ""
echo "--- Preparing BLAST database for ${SPECIES} proteome ---"
if [[ ! -f "${QUERY_PROTEOME}.pdb" ]]; then
    makeblastdb \
        -in "${QUERY_PROTEOME}" \
        -dbtype prot \
        -out "${QUERY_PROTEOME}" \
        -parse_seqids \
        -title "${SPECIES}_proteome"
    echo "  BLAST database created."
else
    echo "  BLAST database already exists. Skipping."
fi

# ── Generate gene-to-protein CSV ─────────────────────────────
echo ""
echo "--- Generating gene-to-protein mapping ---"
GENE_TO_PROTEIN_CSV="${RUN_DIR}/${SPECIES}_gene_to_protein.csv"
python3 << 'PYEOF'
import re, csv
fasta = "${QUERY_PROTEOME}"
out   = "${GENE_TO_PROTEIN_CSV}"
rows  = []
with open(fasta) as f:
    for line in f:
        if not line.startswith('>'):
            continue
        line = line[1:].strip()
        parts = line.split()
        gene_id = parts[0]
        desc = ' '.join(parts[1:])
        # Extract locus tag from NCBI headers (e.g. P5673_033792)
        # For non-NCBI species, locus_tag = gene_id (self-mapping)
        match = re.search(r'\b([A-Z][A-Z0-9]+_\d{4,})\b', desc)
        locus_tag = match.group(1) if match else gene_id
        rows.append((gene_id, locus_tag))
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene_id', 'protein_id'])
    w.writerows(rows)
print(f"  Written: {out} ({len(rows)} entries)")
PYEOF
# Only use auto-generated CSV if no pre-configured path exists
[[ -z "${GENE_TO_PROTEIN}" ]] && export GENE_TO_PROTEIN="${GENE_TO_PROTEIN_CSV}"

# ── Assemble OrthoFinder input directory ────────────────────
echo ""
echo "--- Assembling OrthoFinder input directory ---"
mkdir -p "${OF_INPUT_DIR}"
ln -sf "${QUERY_PROTEOME}" "${OF_INPUT_DIR}/${SPECIES}.faa" 2>/dev/null || true
ln -sf "${HUMAN_FASTA}"    "${OF_INPUT_DIR}/Hsap.faa"       2>/dev/null || true
ln -sf "${MOUSE_FASTA}"    "${OF_INPUT_DIR}/Mmus.faa"       2>/dev/null || true
ln -sf "${NVEC_PROTEOME}"  "${OF_INPUT_DIR}/Nvec.faa"       2>/dev/null || true

echo "  Input proteomes:"
for f in "${OF_INPUT_DIR}"/*.faa; do
    echo "    $(basename ${f}): $(grep -c '^>' ${f}) sequences"
done

# ── Chunk proteome for DeepTMHMM ────────────────────────────
echo ""
echo "--- Chunking proteome for DeepTMHMM (2000 seqs/chunk) ---"
python3 << PYEOF
from pathlib import Path
fasta = "${QUERY_PROTEOME}"
out_dir = Path("${RUN_DIR}/05_tmhmm/chunks")
out_dir.mkdir(parents=True, exist_ok=True)
records = []
current_header, current_seq = None, []
with open(fasta) as f:
    for line in f:
        if line.startswith('>'):
            if current_header:
                records.append((current_header, ''.join(current_seq)))
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
if current_header:
    records.append((current_header, ''.join(current_seq)))
chunk_size = 2000
chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]
print(f"  Total sequences: {len(records):,}")
print(f"  Chunks ({chunk_size} seqs each): {len(chunks)}")
with open("${RUN_DIR}/05_tmhmm/n_chunks.txt", "w") as f:
    f.write(str(len(chunks)))
for i, chunk in enumerate(chunks):
    out_file = out_dir / f"chunk_{i+1:02d}.fasta"
    if not out_file.exists():
        with open(out_file, 'w') as f:
            for header, seq in chunk:
                f.write(header)
                f.write(''.join(seq))
        print(f"  Written: {out_file.name} ({len(chunk)} seqs)")
    else:
        print(f"  EXISTS:  {out_file.name}")
PYEOF

N_CHUNKS=$(cat "${RUN_DIR}/05_tmhmm/n_chunks.txt")
echo "  DeepTMHMM chunks: ${N_CHUNKS}"

# ── Helper: submit a BSUB script and return numeric JID only ─
submit_job() {
    local job_name="$1"
    local script_file="$2"
    local raw_output jid
    raw_output=$(bsub < "${script_file}" 2>&1)
    jid=$(echo "${raw_output}" | grep -oP '(?<=Job <)\d+' || true)
    echo "  Submitted ${job_name}: JID=${jid}" >&2
    echo "${jid}"
}

echo ""
echo "--- Generating and submitting job scripts ---"

# ============================================================
# JOB 01: eggNOG-mapper
# ============================================================
cat > "${LOG_DIR}/job_01_eggnog.bsub" << BSUB
#!/bin/bash
#BSUB -J eggnog_${SPECIES}
#BSUB -q bigmem
#BSUB -P ${LSF_PROJECT}
#BSUB -n 16
#BSUB -W 12:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o ${LOG_DIR}/01_eggnog_%J.out
#BSUB -e ${LOG_DIR}/01_eggnog_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N
set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env
export EGGNOG_DATA_DIR="${EGGNOG_DB}"
echo "=== 01_eggnog: ${SPECIES} ==="
date && hostname
mkdir -p "${EGGNOG_OUT}"
emapper.py \
    -i "${QUERY_PROTEOME}" \
    --itype proteins \
    -o "${SPECIES}" \
    --output_dir "${EGGNOG_OUT}" \
    --data_dir "${EGGNOG_DB}" \
    --tax_scope Metazoa \
    --go_evidence non-electronic \
    --target_orthologs all \
    --cpu 16 \
    --override
echo "=== eggNOG complete ===" && date
BSUB

# ============================================================
# JOB 02: OrthoFinder
# ============================================================
cat > "${LOG_DIR}/job_02_orthofinder.bsub" << BSUB
#!/bin/bash
#BSUB -J orthofinder_${SPECIES}
#BSUB -q bigmem
#BSUB -P ${LSF_PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o ${LOG_DIR}/02_orthofinder_%J.out
#BSUB -e ${LOG_DIR}/02_orthofinder_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate orthofinder_env

echo "=== 02_orthofinder: ${SPECIES} ==="
date && hostname

orthofinder \
    -f "${OF_INPUT_DIR}" \
    -o "${OF_OUT}/Results_OrthoFinder" \
    -t 16 -a 16 -M msa -A mafft -T fasttree

echo "=== OrthoFinder complete ===" && date
BSUB

# ============================================================
# JOB 03: InterProScan
# ============================================================
cat > "${LOG_DIR}/job_03_interproscan.bsub" << BSUB
#!/bin/bash
#BSUB -J ips_${SPECIES}
#BSUB -q bigmem
#BSUB -P ${LSF_PROJECT}
#BSUB -n 16
#BSUB -W 48:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o ${LOG_DIR}/03_interproscan_%J.out
#BSUB -e ${LOG_DIR}/03_interproscan_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env

echo "=== 03_interproscan: ${SPECIES} ==="
date && hostname
mkdir -p "${IPS_OUT}/tmp"

${INTERPROSCAN_BIN} \
    -i "${QUERY_PROTEOME}" \
    -f tsv -goterms -pa \
    -o "${IPS_TSV}" \
    --cpu 16 \
    --tempdir "${IPS_OUT}/tmp"

echo "=== InterProScan complete ===" && date
BSUB

# ============================================================
# JOB 04: RBH vs SwissProt
# ============================================================
cat > "${LOG_DIR}/job_04_rbh.bsub" << BSUB
#!/bin/bash
#BSUB -J rbh_${SPECIES}
#BSUB -q bigmem
#BSUB -P ${LSF_PROJECT}
#BSUB -n 16
#BSUB -W 08:00
#BSUB -R "rusage[mem=4000]"
#BSUB -o ${LOG_DIR}/04_rbh_%J.out
#BSUB -e ${LOG_DIR}/04_rbh_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env

echo "=== 04_rbh: ${SPECIES} ==="
date && hostname
mkdir -p "${RBH_OUT}"

blastp -query "${QUERY_PROTEOME}" -db "${SWISSPROT_FASTA}" \
    -out "${RBH_OUT}/${SPECIES}_vs_sprot.blast6" \
    -outfmt "6 qseqid sseqid pident evalue bitscore" \
    -evalue 1e-5 -max_target_seqs 1 -num_threads 16

blastp -query "${SWISSPROT_FASTA}" -db "${QUERY_PROTEOME}" \
    -out "${RBH_OUT}/sprot_vs_${SPECIES}.blast6" \
    -outfmt "6 qseqid sseqid pident evalue bitscore" \
    -evalue 1e-5 -max_target_seqs 1 -num_threads 16

python3 "${RUN_DIR}/py_helpers/parse_rbh.py" \
    "${SPECIES}" "${RBH_OUT}" "${RBH_TSV}" "${SWISSPROT_FASTA}"

echo "=== RBH complete ===" && date
BSUB

# ============================================================
# JOB 05: SignalP 6.0
# ============================================================
cat > "${LOG_DIR}/job_05_signalp.bsub" << BSUB
#!/bin/bash
#BSUB -J signalp_${SPECIES}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 8
#BSUB -W 20:00
#BSUB -R "rusage[mem=16000]"
#BSUB -o ${LOG_DIR}/05_signalp_%J.out
#BSUB -e ${LOG_DIR}/05_signalp_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env
export LD_LIBRARY_PATH="/nethome/kxw755/anaconda3/envs/annotation_env/lib:\${LD_LIBRARY_PATH}"

echo "=== 05_signalp: ${SPECIES} ==="
date && hostname
mkdir -p "${SIGNALP_OUT}"

if ! command -v signalp6 &>/dev/null; then
    echo "WARNING: signalp6 not found. Skipping."
    echo "SKIPPED" > "${SIGNALP_OUT}/signalp_not_installed.txt"
    exit 0
fi

# Strip long FASTA headers to avoid OSError in SignalP plot generation
SIGNALP_FASTA="\${SIGNALP_OUT}/input_clean.fasta"
python3 -c "
with open('${QUERY_PROTEOME}') as fi, open('\${SIGNALP_FASTA}', 'w') as fo:
    for line in fi:
        fo.write(('>' + line[1:].split()[0] + '\n') if line.startswith('>') else line)
print('Clean FASTA:', '\${SIGNALP_FASTA}')
"

signalp6 \
    --fastafile "\${SIGNALP_FASTA}" \
    --output_dir "${SIGNALP_OUT}" \
    --format txt \
    --organism eukarya \
    --mode slow-sequential \
    --torch_num_threads 8

echo "=== SignalP complete ===" && date
BSUB

# ============================================================
# JOB 05b: DeepTMHMM (parallel chunks)
# ============================================================
TMHMM_DIR="${RUN_DIR}/05_tmhmm"
DTMHMM_DIR="/nethome/kxw755/DeepTMHMM-Academic-License-v1.0"
CHUNK_JIDS=()

for i in $(seq -w 1 "${N_CHUNKS}"); do
    CHUNK="${TMHMM_DIR}/chunks/chunk_${i}.fasta"
    CHUNK_OUT="${TMHMM_DIR}/results/chunk_${i}"

    cat > "${LOG_DIR}/job_05b_tmhmm_${i}.bsub" << BSUB
#!/bin/bash
#BSUB -J tmhmm_${SPECIES}_${i}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 4
#BSUB -W 12:00
#BSUB -R "rusage[mem=16000]"
#BSUB -o ${LOG_DIR}/05b_tmhmm_${i}_%J.out
#BSUB -e ${LOG_DIR}/05b_tmhmm_${i}_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env
export LD_LIBRARY_PATH="/nethome/kxw755/anaconda3/envs/annotation_env/lib:\${LD_LIBRARY_PATH}"

echo "=== DeepTMHMM chunk ${i}: ${SPECIES} ==="
date && hostname

rm -rf "${CHUNK_OUT}"
cd "${DTMHMM_DIR}"
python3 predict.py --fasta "${CHUNK}" --output "${CHUNK_OUT}"

echo "=== DeepTMHMM chunk ${i} complete ===" && date
ls -lh "${CHUNK_OUT}/"
BSUB

    JID=$(submit_job "05b_tmhmm_${i}" "${LOG_DIR}/job_05b_tmhmm_${i}.bsub")
    CHUNK_JIDS+=("${JID}")
done

# ── Concatenate DeepTMHMM chunks after all complete ──────────
TMHMM_DEP=$(printf "done(%s) && " "${CHUNK_JIDS[@]}")
TMHMM_DEP="${TMHMM_DEP% && }"

cat > "${LOG_DIR}/job_05b_tmhmm_merge.bsub" << BSUB
#!/bin/bash
#BSUB -J tmhmm_merge_${SPECIES}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 1
#BSUB -W 00:30
#BSUB -R "rusage[mem=4000]"
#BSUB -w "${TMHMM_DEP}"
#BSUB -o ${LOG_DIR}/05b_tmhmm_merge_%J.out
#BSUB -e ${LOG_DIR}/05b_tmhmm_merge_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
echo "=== Merging DeepTMHMM chunks: ${SPECIES} ==="
date

TMHMM_DIR="${TMHMM_DIR}"
OUT_GFF3="\${TMHMM_DIR}/TMRs.gff3"

echo "##gff-version 3" > "\${OUT_GFF3}"
for chunk_dir in "${TMHMM_DIR}/results"/chunk_*/; do
    if [[ -f "\${chunk_dir}/TMRs.gff3" ]]; then
        grep -v "^##gff-version" "\${chunk_dir}/TMRs.gff3" >> "\${OUT_GFF3}"
    fi
done

N_GENES=\$(grep -c "^# " "\${OUT_GFF3}" || true)
echo "Merged: \${N_GENES} sequences in \${OUT_GFF3}"
echo "=== DeepTMHMM merge complete ===" && date
BSUB

JID_05b_MERGE=$(submit_job "05b_tmhmm_merge" "${LOG_DIR}/job_05b_tmhmm_merge.bsub")

# ============================================================
# JOB 06: Tier 3 HMM search
# ============================================================
cat > "${LOG_DIR}/job_06_tier3.bsub" << BSUB
#!/bin/bash
#BSUB -J tier3_${SPECIES}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 8
#BSUB -W 06:00
#BSUB -R "rusage[mem=8000]"
#BSUB -o ${LOG_DIR}/06_tier3_%J.out
#BSUB -e ${LOG_DIR}/06_tier3_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env

echo "=== 06_tier3: ${SPECIES} ==="
date && hostname
mkdir -p "${TIER3_OUT}"

hmmsearch \
    --cpu 8 \
    --domtblout "${TIER3_DOMTBL}" \
    -E 1e-3 --domE 1e-3 \
    "${TIER3_HMM}" \
    "${QUERY_PROTEOME}" \
    > "${TIER3_OUT}/${SPECIES}_tier3_hmmsearch.log"

echo "Hits: \$(grep -v '^#' ${TIER3_DOMTBL} | wc -l)"
echo "=== Tier 3 complete ===" && date
BSUB

# ============================================================
# JOB 07: Tier 4 BLASTp
# ============================================================
cat > "${LOG_DIR}/job_07_tier4.bsub" << BSUB
#!/bin/bash
#BSUB -J tier4_${SPECIES}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 8
#BSUB -W 01:00
#BSUB -R "rusage[mem=4000]"
#BSUB -o ${LOG_DIR}/07_tier4_%J.out
#BSUB -e ${LOG_DIR}/07_tier4_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env

echo "=== 07_tier4: ${SPECIES} ==="
date && hostname
mkdir -p "${TIER4_OUT}"

blastp \
    -query  "${QUERY_PROTEOME}" \
    -db     "${TIER4_FASTA}" \
    -out    "${TIER4_OUT}/${SPECIES}_vs_tier4_raw.blast6" \
    -outfmt "6 qseqid sseqid pident evalue bitscore" \
    -evalue 1e-5 -max_target_seqs 5 -num_threads 8

python3 "${RUN_DIR}/py_helpers/parse_tier4.py" \
    "${SPECIES}" "${TIER4_OUT}" "${TIER4_BLAST}"

echo "=== Tier 4 complete ===" && date
BSUB

# ── Write Python helper scripts ──────────────────────────────
cat > "${RUN_DIR}/py_helpers/parse_rbh.py" << 'PYEOF'
#!/usr/bin/env python3
"""Parse reciprocal best BLAST hits against SwissProt."""
import sys, csv
from pathlib import Path

species, rbh_dir, rbh_tsv, sprot_fasta = sys.argv[1], Path(sys.argv[2]), sys.argv[3], sys.argv[4]
gene_id_col = f"{species.lower()}_gene_id"

fwd_file = rbh_dir / f"{species}_vs_sprot.blast6"
rev_file = rbh_dir / f"sprot_vs_{species}.blast6"

fwd = {}
with open(fwd_file) as f:
    for row in csv.reader(f, delimiter='\t'):
        q, s, pct, ev, score = row
        if q not in fwd:
            fwd[q] = (s, float(pct), float(ev))

rev = {}
with open(rev_file) as f:
    for row in csv.reader(f, delimiter='\t'):
        q, s, pct, ev, score = row
        # Strip database prefix from subject ID (e.g. gb|KAK2546616.1| -> KAK2546616.1)
        s_clean = s.split('|')[1] if s.count('|') >= 2 else s
        if q not in rev:
            rev[q] = (s_clean, float(pct), float(ev))

sprot_names = {}
with open(sprot_fasta) as f:
    for line in f:
        if not line.startswith('>'):
            continue
        parts = line[1:].split()
        if '|' in parts[0]:
            tokens = parts[0].split('|')
            if len(tokens) >= 3:
                sprot_names[tokens[1]] = tokens[2]
                sprot_names[tokens[2]] = tokens[2]

rbh_rows = []
for qgene, (sprot_hit, pct, evalue) in fwd.items():
    if '|' in sprot_hit:
        tokens = sprot_hit.split('|')
        sprot_acc   = tokens[1] if len(tokens) > 1 else sprot_hit
        sprot_entry = tokens[2] if len(tokens) > 2 else sprot_hit
    else:
        sprot_acc   = sprot_hit
        sprot_entry = sprot_names.get(sprot_hit, sprot_hit)
    rev_best = rev.get(sprot_hit, (None,))[0] or rev.get(sprot_acc, (None,))[0]
    if rev_best == qgene:
        rbh_rows.append({
            gene_id_col:         qgene,
            'rbh_swissprot_acc': sprot_acc,
            'rbh_gene_name':     sprot_names.get(sprot_acc, sprot_entry),
            'rbh_pct_identity':  round(pct, 2),
            'rbh_evalue':        evalue,
        })

print(f"Forward: {len(fwd)} | Reverse: {len(rev)} | RBH: {len(rbh_rows)}")
with open(rbh_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t',
        fieldnames=[gene_id_col, 'rbh_swissprot_acc', 'rbh_gene_name',
                    'rbh_pct_identity', 'rbh_evalue'])
    w.writeheader()
    w.writerows(rbh_rows)
print(f"Written: {rbh_tsv}")
PYEOF

cat > "${RUN_DIR}/py_helpers/parse_tier4.py" << 'PYEOF'
#!/usr/bin/env python3
"""Parse Tier 4 BLAST hits."""
import sys, csv
from pathlib import Path

species, tier4_dir, tier4_out = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
gene_id_col = f"{species.lower()}_gene_id"
raw_file = tier4_dir / f"{species}_vs_tier4_raw.blast6"

rows = []
with open(raw_file) as f:
    for line in csv.reader(f, delimiter='\t'):
        qseqid, sseqid, pident, evalue, bitscore = line
        group = sseqid.split('_')[0] if '_' in sseqid else 'Unknown'
        rows.append({gene_id_col: qseqid, 'tier4_blast_hit': sseqid,
                     'tier4_source': group, 'tier4_pct_identity': float(pident),
                     'tier4_evalue': float(evalue), 'tier4_bitscore': float(bitscore)})

seen, dedup = set(), []
for row in sorted(rows, key=lambda r: r['tier4_evalue']):
    key = (row[gene_id_col], row['tier4_source'])
    if key not in seen:
        dedup.append(row)
        seen.add(key)

with open(tier4_out, 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t',
        fieldnames=[gene_id_col, 'tier4_blast_hit', 'tier4_source',
                    'tier4_pct_identity', 'tier4_evalue', 'tier4_bitscore'])
    w.writeheader()
    w.writerows(dedup)
print(f"Tier 4 hits: {len(dedup)}")
print(f"Written: {tier4_out}")
PYEOF

# ── Submit jobs 01–07 ────────────────────────────────────────
JID_01=$(submit_job "01_eggnog"       "${LOG_DIR}/job_01_eggnog.bsub")
JID_02=$(submit_job "02_orthofinder"  "${LOG_DIR}/job_02_orthofinder.bsub")
JID_03=$(submit_job "03_interproscan" "${LOG_DIR}/job_03_interproscan.bsub")
JID_04=$(submit_job "04_rbh"          "${LOG_DIR}/job_04_rbh.bsub")
JID_05=$(submit_job "05_signalp"      "${LOG_DIR}/job_05_signalp.bsub")
JID_06=$(submit_job "06_tier3"        "${LOG_DIR}/job_06_tier3.bsub")
JID_07=$(submit_job "07_tier4"        "${LOG_DIR}/job_07_tier4.bsub")

# Save JIDs
printf "%s\n" "${JID_01}" "${JID_02}" "${JID_03}" "${JID_04}" \
              "${JID_05}" "${JID_05b_MERGE}" "${JID_06}" "${JID_07}" \
    > "${LOG_DIR}/submitted_jids.txt"

# ── Job 08: Merge — depends on all jobs including DeepTMHMM ─
echo ""
echo "--- Building merge job ---"

DEPENDENCY="done(${JID_01}) && done(${JID_02}) && done(${JID_03}) && \
done(${JID_04}) && done(${JID_05}) && done(${JID_05b_MERGE}) && \
done(${JID_06}) && done(${JID_07})"

cat > "${LOG_DIR}/job_08_merge.bsub" << BSUB
#!/bin/bash
#BSUB -J merge_${SPECIES}
#BSUB -q general
#BSUB -P ${LSF_PROJECT}
#BSUB -n 4
#BSUB -W 02:00
#BSUB -R "rusage[mem=8000]"
#BSUB -w "${DEPENDENCY}"
#BSUB -o ${LOG_DIR}/08_merge_%J.out
#BSUB -e ${LOG_DIR}/08_merge_%J.err
#BSUB -u ${EMAIL}
#BSUB -B
#BSUB -N

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate annotation_env
export LD_LIBRARY_PATH="/nethome/kxw755/anaconda3/envs/annotation_env/lib:\${LD_LIBRARY_PATH}"

echo "=== 08_merge: ${SPECIES} ==="
date && hostname
mkdir -p "${FINAL_OUT}/gene_lists"

python3 "/scratch/dark_genes/annotation_pipeline/scripts/08_merge_annotate.py" \
    --species         "${SPECIES}" \
    --run_dir         "${RUN_DIR}" \
    --proteome_fasta  "${QUERY_PROTEOME}" \
    --gff3            "${SPECIES_GFF3}" \
    --eggnog          "${EGGNOG_ANNOT}" \
    --orthofinder     "${OF_RESULTS_DIR}" \
    --interproscan    "${IPS_TSV}" \
    --rbh             "${RBH_TSV}" \
    --signalp         "${SIGNALP_TSV}" \
    --tmhmm           "${RUN_DIR}/05_tmhmm/TMRs.gff3" \
    --tier3           "${TIER3_DOMTBL}" \
    --tier4           "${TIER4_BLAST}" \
    --out_table       "${MASTER_TABLE}" \
    --out_genelists   "${GENE_LISTS_DIR}" \
    --genome_version  "${GENOME_VERSION}" \
    --proteome_source "${PROTEOME_SOURCE}" \
    --staging_dir     "${STAGING_DIR}" \
    --gene_to_protein "${GENE_TO_PROTEIN}"

echo "=== Pipeline complete for ${SPECIES} ===" && date
BSUB

JID_08=$(submit_job "08_merge" "${LOG_DIR}/job_08_merge.bsub")
echo "${JID_08}" >> "${LOG_DIR}/submitted_jids.txt"

# ── Summary ──────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " All jobs submitted for species: ${SPECIES}"
echo "============================================================"
echo "  01_eggnog        : JID ${JID_01}"
echo "  02_orthofinder   : JID ${JID_02}"
echo "  03_interproscan  : JID ${JID_03}"
echo "  04_rbh           : JID ${JID_04}"
echo "  05_signalp       : JID ${JID_05}  (~8-20hr)"
echo "  05b_tmhmm chunks : JIDs ${CHUNK_JIDS[*]}"
echo "  05b_tmhmm merge  : JID ${JID_05b_MERGE}  (waits on chunks)"
echo "  06_tier3         : JID ${JID_06}"
echo "  07_tier4         : JID ${JID_07}"
echo "  08_merge         : JID ${JID_08}  (waits on all above)"
echo ""
echo "  Run directory    : ${RUN_DIR}"
echo "  Monitor          : bjobs -u kxw755 | grep ${SPECIES}"
echo "  Logs             : ${LOG_DIR}/"
echo "  Final output     : ${MASTER_TABLE}"
echo "============================================================"
date
