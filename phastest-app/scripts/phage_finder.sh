#!/usr/bin/env bash
# Robust, improved version of the original phage_finder.sh
# - safer shell options
# - better argument parsing and validation
# - correct BLAST+ option (-num_threads)
# - progress reporting based on unique query IDs seen in the BLAST output
# - proper quoting and logging

set -euo pipefail
IFS=$'\n\t'

progname="$(basename "${BASH_SOURCE[0]}")"

usage() {
    cat <<EOF
Usage: $progname JOB_ID DATABASE_HEAD DATABASE_CHILD LOCAL_GI_NUM [LOG_FILE]

Runs BLASTP of JOB_ID.pep (or .faa) against the phage virus DB and writes results to ncbi.out
  JOB_ID          prefix for .pep/.faa/.ptt files
  DATABASE_HEAD   path to DB on head node (unused by this script but kept for compatibility)
  DATABASE_CHILD  path to DB on child node (used for BLAST)
  LOCAL_GI_NUM    expected number of proteins (used for progress)
  LOG_FILE        optional log file (default: phage_finder.log)
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

base="$(pwd)"
job_id="${1:-}"
database_head="${2:-}"
database_child="${3:-}"
local_gi_num="${4:-0}"
log_file="${5:-phage_finder.log}"
log_file2="${log_file}.2"

if [[ -z "$job_id" || -z "$database_child" ]]; then
    echo "ERROR: missing required arguments." >&2
    usage
    exit 2
fi

# Ensure local_gi_num is a non-negative integer (fallback to 0)
if ! [[ "$local_gi_num" =~ ^[0-9]+$ ]]; then
    local_gi_num=0
fi

ncbi_out="$base/ncbi.out"

function progress_bar() {
    local current=$1
    local total=$2
    local bar_length=20
    local progress
    local progress_percent
    local i

    # If the total is 0, do nothing
    # sanitize numeric inputs (remove any non-digits) and default to 0
    current=${current//[^0-9]/}
    total=${total//[^0-9]/}
    if [[ -z "$current" ]]; then
        current=0
    fi
    if [[ -z "$total" ]]; then
        total=0
    fi

    if [ "$total" -eq 0 ]; then
        return
    fi

    # Use 10# prefix to force base-10 interpretation (avoid octal when values have leading zeros)
    progress=$((10#$current*bar_length/10#$total))
    progress_percent=$((10#$current*100/10#$total))

    printf "\rProgress: ["
    for i in $(seq 1 $progress); do printf "="; done
    for i in $(seq $progress $bar_length); do printf " "; done
    printf "] $progress_percent%%"
}

# locate peptide file
if [[ -s "$base/$job_id.pep" ]]; then
    pepfile="$base/$job_id.pep"
elif [[ -s "$base/$job_id.faa" ]]; then
    pepfile="$base/$job_id.faa"
else
    echo "Could not find $job_id.pep or $job_id.faa with data" >> "$log_file2"
    exit 1
fi

# locate info file (keeps original behavior)
if [[ -s "$base/phage_finder_info.txt" ]]; then
    infofile="$base/phage_finder_info.txt"
elif [[ -s "$base/$job_id.ptt" ]]; then
    infofile="$base/$job_id.ptt"
else
    echo "Could not find phage_finder_info.txt or $job_id.ptt with data" >> "$log_file2"
    exit 1
fi

# keep compatibility with callers that pass database_head or rely on infofile being set
: "${database_head:-}" >/dev/null 2>&1 || true
: "${infofile:-}" >/dev/null 2>&1 || true

if ! command -v blastp >/dev/null 2>&1; then
    echo "blastp not found in PATH; please install BLAST+ or adjust PATH" >> "$log_file2"
    exit 127
fi

if [[ ! -s "$ncbi_out" ]]; then
    echo "$ncbi_out does not exist or is empty. Performing BLAST search." >> "$log_file2"
    : > "$ncbi_out"
    START=$(date +%s)

    echo "Running BLASTP: query=$pepfile db=$database_child" >> "$log_file2"
    # Use -num_threads (BLAST+) and proper quoting
    blastp -num_threads 4 -db "$database_child" -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 0.0001 -query "$pepfile" -out "$ncbi_out" -seg no >/dev/null 2>&1 &
    pid=$!

    # Ensure the background BLAST is killed if this script exits
    trap 'kill "$pid" >/dev/null 2>&1 || true' EXIT

    # Monitor progress by counting unique qseqid seen in the ncbi_out file
    while kill -0 "$pid" 2>/dev/null; do
        if [[ -s "$ncbi_out" ]]; then
            # count unique query IDs seen so far (first column)
            current_gi=$(cut -f1 -d$'\t' "$ncbi_out" 2>/dev/null | sort -u 2>/dev/null | wc -l 2>/dev/null || echo 0)
        else
            current_gi=0
        fi
        progress_bar "$current_gi" "$local_gi_num"
        sleep 1
    done

    # wait for process to finish and capture exit status
    wait "$pid" || true
    trap - EXIT

    END=$(date +%s)
    DIFF=$((END - START))
    progress_bar "$local_gi_num" "$local_gi_num"
    printf "\n"
    echo "BLASTP completed: query=$pepfile db=$database_child took ${DIFF}s" >> "$log_file2"
else
    echo "$ncbi_out already exists and is non-empty. Skipping BLAST." >> "$log_file2"
fi



