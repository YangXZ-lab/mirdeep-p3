#!/usr/bin/env bash
# =============================================================================
# mirdeep-p3 test.sh — End-to-end smoke test for source installs
#
# Runs the same identification + annotation pipeline used by GitHub Actions CI
# on a small Arabidopsis test dataset (tests/data/mini.fq + mini_genome.fasta).
#
# Usage:
#   ./test.sh            # run full test (identification + annotation)
#   ./test.sh -q         # quiet mode (suppress pipeline stdout)
#   ./test.sh -k         # keep test output (default: removed on success)
#
# Exit codes:
#   0  all tests passed
#   1  dependency check failed
#   2  data download failed
#   3  identification failed
#   4  annotation failed
# =============================================================================
set -u

# ---- Configuration ----------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE_TAG="mirdeep-p3-v3.1.4c-full"
INDEX_URL="https://github.com/YangXZ-lab/mirdeep-p3/releases/download/${RELEASE_TAG}/data-index.tar.gz"
INPUT_FASTQ="${ROOT}/tests/data/mini.fq"
GENOME="${ROOT}/tests/data/mini_genome.fasta"
OUT_DIR="${ROOT}/tests/test_output"
THREADS=2
KEEP_OUTPUT=0
QUIET=0

# ---- Parse args -------------------------------------------------------------
while getopts "qkh" opt; do
    case "$opt" in
        q) QUIET=1 ;;
        k) KEEP_OUTPUT=1 ;;
        h)
            grep "^#" "$0" | sed 's/^# \{0,1\}//' | head -25
            exit 0
            ;;
        *) exit 1 ;;
    esac
done

# ---- Helpers ----------------------------------------------------------------
log()  { echo "[test] $*"; }
fail() { echo "[test] ERROR: $*" >&2; exit "${2:-1}"; }

run_quiet() {
    if [ "$QUIET" -eq 1 ]; then
        "$@" > /dev/null 2>&1
    else
        "$@"
    fi
}

# ---- 1. Dependency check ----------------------------------------------------
log "Checking dependencies..."
MISSING=()
for tool in bowtie RNAfold samtools bedtools seqkit blastn; do
    command -v "$tool" >/dev/null 2>&1 || MISSING+=("$tool")
done
if [ "${#MISSING[@]}" -gt 0 ]; then
    fail "Missing external tools: ${MISSING[*]} (activate your conda env first)" 1
fi
log "All required tools found."

# ---- 2. Test data -----------------------------------------------------------
if [ ! -f "$INPUT_FASTQ" ] || [ ! -f "$GENOME" ]; then
    fail "Test data not found in ${ROOT}/tests/data/" 1
fi

# ---- 3. PmiREN / Rfam index (needed by annotation & identification) ---------
if [ ! -d "${ROOT}/data/index" ]; then
    log "data/index not found - downloading from GitHub Release..."
    mkdir -p "${ROOT}/data"
    if ! wget -q -O /tmp/data-index.tar.gz "$INDEX_URL" \
        && ! curl -sL -o /tmp/data-index.tar.gz "$INDEX_URL"; then
        rm -f /tmp/data-index.tar.gz
        fail "Could not download ${INDEX_URL}" 2
    fi
    if ! tar xzf /tmp/data-index.tar.gz -C "${ROOT}/data/"; then
        rm -f /tmp/data-index.tar.gz
        fail "Could not extract data-index.tar.gz" 2
    fi
    rm -f /tmp/data-index.tar.gz
    log "Index data extracted to ${ROOT}/data/index/"
else
    log "Using existing ${ROOT}/data/index/"
fi

# ---- 4. Clean previous output -----------------------------------------------
rm -rf "$OUT_DIR"

# ---- 5. Identification ------------------------------------------------------
log "Running identification (this may take a few minutes)..."
if ! run_quiet "${ROOT}/mirdeep-p3" identification \
        -i "$INPUT_FASTQ" \
        -o "$OUT_DIR" \
        -g "$GENOME" \
        -t "$THREADS"; then
    fail "Identification step failed (see log above)" 3
fi

# ---- 6. Verify identification outputs ---------------------------------------
ID_OUT="$OUT_DIR/mini"
for f in "$ID_OUT/mini-processed.fa" "$ID_OUT/mini_filter_P_prediction"; do
    [ -s "$f" ] || fail "Identification output missing: $f" 3
done
log "Identification outputs verified."

# ---- 7. Annotation ----------------------------------------------------------
log "Running annotation..."
if ! run_quiet "${ROOT}/mirdeep-p3" annotation \
        -i "$ID_OUT" \
        -g "$GENOME" \
        -d "$ID_OUT/index/genome_index" \
        -t "$THREADS" \
        --prefix_miRNA "Ath" \
        --prefix mirdp3-test \
        --species "Arabidopsis thaliana" \
        -o "$OUT_DIR/annotation"; then
    fail "Annotation step failed (see log above)" 4
fi

# ---- 8. Verify annotation outputs -------------------------------------------
ANNO_OUT="$OUT_DIR/annotation/mirdp3-test"
if [ ! -s "$ANNO_OUT/mirdp3-test-basic-info" ]; then
    fail "Annotation output missing: $ANNO_OUT/mirdp3-test-basic-info" 4
fi
log "Annotation outputs verified."

# ---- 9. Cleanup -------------------------------------------------------------
if [ "$KEEP_OUTPUT" -eq 0 ]; then
    log "Removing test output (use -k to keep)..."
    rm -rf "$OUT_DIR"
fi

# ---- Done -------------------------------------------------------------------
log "=============================================="
log " ALL TESTS PASSED "
log "=============================================="
exit 0
