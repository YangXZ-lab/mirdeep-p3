#!/usr/bin/env bash
# =============================================================================
# mirdeep-p3 test.sh — End-to-end smoke test for source installs
#
# Runs the same identification + annotation pipeline used by GitHub Actions CI
# on a small Arabidopsis test dataset (tests/data/mini.fq + mini_genome.fasta).
#
# Usage:
#   ./test.sh            # run full test (identification + annotation + target prediction
#                        + analysis smoke)
#   ./test.sh -q         # quiet mode (suppress pipeline stdout)
#   ./test.sh -k         # keep test output (default: removed on success)
#   ./test.sh -t 16      # use 16 threads (default: 8)
#
# Exit codes:
#   0  all tests passed
#   1  dependency check failed
#   2  data download failed
#   3  identification failed
#   4  annotation failed
#   5  target prediction failed
#   6  analysis smoke test failed (Stat / Differential_expression)
# =============================================================================
set -u

# ---- Configuration ----------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The Release that hosts data-index.tar.gz is recorded in data-index.env -- a single
# source of truth shared with the GitHub workflows, so the index can be updated
# without touching the software version.
RELEASE_TAG=""
if [ -f "${ROOT}/data-index.env" ]; then
    RELEASE_TAG="$(sed -n 's/^DATA_INDEX_RELEASE=//p' "${ROOT}/data-index.env" | head -1 | tr -d '[:space:]')"
fi
RELEASE_TAG="${RELEASE_TAG:-mirdeep-p3-v3.1.4c-full}"
INDEX_URL="https://github.com/YangXZ-lab/mirdeep-p3/releases/download/${RELEASE_TAG}/data-index.tar.gz"
INPUT_FASTQ="${ROOT}/tests/data/mini.fq"
GENOME="${ROOT}/tests/data/mini_genome.fasta"
OUT_DIR="${ROOT}/tests/test_output"
THREADS=8
KEEP_OUTPUT=0
QUIET=0

# ---- Parse args -------------------------------------------------------------
while getopts "t:qkh" opt; do
    case "$opt" in
        t) THREADS="$OPTARG" ;;
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
    # Record the release that was just unpacked so the software can report it (O-10).
    printf '%s\n' "$RELEASE_TAG" > "${ROOT}/data/index/.index-version" 2>/dev/null || true
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

# ---- 9. Target prediction ---------------------------------------------------
# Fixture: 20 mature miRNAs and 300 CDS extracted from the Arctium_lappa dataset
# used to validate the psrna engine (docs/PSRNA-ENGINE.md). It is small enough to
# run in CI and fixed enough to assert an exact row count.
TF_MATURE="$ROOT/tests/data/target_mini_mature.fa"
TF_CDS="$ROOT/tests/data/target_mini_cds.fa"
EXPECTED_PSRNA_ROWS=120
if [ -f "$TF_MATURE" ] && [ -f "$TF_CDS" ]; then
    log "Running target prediction (default engine: psrna)..."
    TF_OUT="$OUT_DIR/target_psrna"
    if ! run_quiet "$ROOT/mirdeep-p3" analysis Target_finder \
            -i "$TF_MATURE" \
            -c "$TF_CDS" \
            -o "$TF_OUT" \
            -t "$THREADS" \
            -e 2.5; then
        fail "Target_finder (psrna) failed (see log above)" 5
    fi
    TF_FILE="$TF_OUT/target_finder.tsv"
    [ -s "$TF_FILE" ] || fail "Target prediction output missing: $TF_FILE" 5
    TF_ROWS="$(wc -l < "$TF_FILE")"
    if [ "$TF_ROWS" -ne "$EXPECTED_PSRNA_ROWS" ]; then
        fail "psrna engine: expected $EXPECTED_PSRNA_ROWS rows, got $TF_ROWS" 5
    fi
    awk -F'\t' 'NF != 6 { exit 1 }' "$TF_FILE" || fail "psrna engine: rows are not 6 columns" 5

    # The output order must be fixed: the same input twice has to be identical.
    TF_OUT2="$OUT_DIR/target_psrna_repeat"
    run_quiet "$ROOT/mirdeep-p3" analysis Target_finder \
        -i "$TF_MATURE" -c "$TF_CDS" -o "$TF_OUT2" -t "$THREADS" -e 2.5 \
        || fail "Target_finder repeat run failed" 5
    cmp -s "$TF_FILE" "$TF_OUT2/target_finder.tsv" \
        || fail "psrna engine: output is not deterministic" 5

    # --no-allow-gap must stay available; on this fixture it changes nothing.
    TF_OUT4="$OUT_DIR/target_psrna_nogap"
    run_quiet "$ROOT/mirdeep-p3" analysis Target_finder --no-allow-gap \
        -i "$TF_MATURE" -c "$TF_CDS" -o "$TF_OUT4" -t "$THREADS" -e 2.5 \
        || fail "Target_finder --no-allow-gap failed" 5
    cmp -s "$TF_FILE" "$TF_OUT4/target_finder.tsv" \
        || fail "psrna engine: --no-allow-gap changed this fixture (expected identical)" 5
    log "Target prediction (psrna) verified: $TF_ROWS rows, deterministic."

    if command -v ssearch36 >/dev/null 2>&1; then
        log "Running target prediction (engine ssearch)..."
        TF_OUT3="$OUT_DIR/target_ssearch"
        if ! run_quiet "$ROOT/mirdeep-p3" analysis Target_finder \
                --engine ssearch \
                -i "$TF_MATURE" \
                -c "$TF_CDS" \
                -o "$TF_OUT3" \
                -t "$THREADS" \
                -e 2.5; then
            fail "Target_finder (ssearch) failed (see log above)" 5
        fi
        TF_ROWS3="$(wc -l < "$TF_OUT3/target_finder.tsv")"
        if [ "$TF_ROWS3" -lt 100 ]; then
            fail "ssearch engine: expected at least 100 rows, got $TF_ROWS3" 5
        fi
        awk -F'\t' 'NF != 6 { exit 1 }' "$TF_OUT3/target_finder.tsv" \
            || fail "ssearch engine: rows are not 6 columns" 5
        log "Target prediction (ssearch) verified: $TF_ROWS3 rows."
    else
        log "ssearch36 not found - skipping the ssearch engine check."
    fi
else
    log "Target prediction fixture not found - skipping target prediction."
fi
# ---- 10. Analysis smoke (Stat + Differential_expression) --------------------
# The two analysis subcommands that need nothing beyond what the annotation step
# already produced.  Skipped rather than failed when Rscript is absent: analysis
# is R-based, and a Python-only environment must not turn the build red.
if command -v Rscript >/dev/null 2>&1; then
    log "Running analysis smoke test (Stat)..."
    STAT_OUT="$OUT_DIR/analysis_stat"
    if ! run_quiet "$ROOT/mirdeep-p3" analysis Stat \
            -i "$ANNO_OUT/mirdp3-test-basic-info" \
            -o "$STAT_OUT"; then
        fail "analysis Stat failed (see log above)" 6
    fi
    for f in length_dist first_base_dist base_dist miRNA_family; do
        [ -s "$STAT_OUT/$f.pdf" ] || fail "Stat output missing: $f.pdf" 6
    done
    log "Stat verified: 4 figures."

    # The DE step rewrites the matrices it is given (filter_low_expression drops
    # low-count miRNAs in place), so it always gets a private copy of the fixture.
    DE_COUNT="$ROOT/tests/data/analysis_mini.count"
    DE_EXP="$ROOT/tests/data/analysis_mini.exp"
    if [ -f "$DE_COUNT" ] && [ -f "$DE_EXP" ]; then
        log "Running analysis smoke test (Differential_expression)..."
        DE_IN="$OUT_DIR/analysis_de_in"
        DE_OUT="$OUT_DIR/analysis_de"
        mkdir -p "$DE_IN"
        cp "$DE_COUNT" "$DE_EXP" "$DE_IN/"
        # Addressed through its alias on purpose: "DE" must reach the same module
        # as the full name (a mis-dispatch used to print the help and exit 1).
        if ! run_quiet "$ROOT/mirdeep-p3" analysis DE \
                -c "$DE_IN/analysis_mini.count" \
                -r "$DE_IN/analysis_mini.exp" \
                --case1 3,4,5 \
                --case2 6,7,8 \
                --case1name mock \
                --case2name treat \
                -o "$DE_OUT"; then
            fail "analysis Differential_expression failed (see log above)" 6
        fi
        DE_RES="$DE_OUT/mock_treat_res.txt"
        for f in "$DE_RES" "$DE_OUT/mock_treat_up.txt" "$DE_OUT/mock_treat_down.txt" \
                 "$DE_OUT/final_miRNA_expression.txt" "$DE_OUT/volcano_mock_treat.pdf"; do
            [ -s "$f" ] || fail "Differential_expression output missing: $f" 6
        done
        DE_ROWS=$(( $(wc -l < "$DE_RES") - 1 ))
        [ "$DE_ROWS" -ge 100 ] \
            || fail "Differential_expression: expected >=100 miRNAs, got $DE_ROWS" 6
        log "Differential expression verified: $DE_ROWS miRNAs tested, DEG lists written."
    else
        log "DE fixture not found - skipping Differential_expression."
    fi
else
    log "Rscript not found - skipping the analysis smoke test."
fi

# ---- 11. CLI, help and unit checks ------------------------------------------
log "Checking generated help, the console entry point and unit tests..."

for combo in "1,3" "2,3" "1,2,3"; do
    HELP_OUT="$("${ROOT}/mirdeep-p3" -h --step "$combo" 2>&1)" || fail "combined help failed for --step $combo" 7
    printf '%s' "$HELP_OUT" | grep -q "Steps: " \
        || fail "combined help for --step $combo has no step list" 7
done
"${ROOT}/mirdeep-p3" -h --step 1,3 2>&1 | grep -q "^analysis subcommands" \
    || fail "--step 1,3 help does not list the analysis subcommands" 7
if "${ROOT}/mirdeep-p3" -h --step 1,3 2>&1 | grep -q -- "--prefix_miRNA"; then
    fail "--step 1,3 help must not list annotation-only options" 7
fi
"${ROOT}/mirdeep-p3" -h --step 2,3 2>&1 | grep -q -- "--prefix_miRNA" \
    || fail "--step 2,3 help should list --prefix_miRNA" 7
"${ROOT}/mirdeep-p3" -h --step 1,2 2>&1 | grep -q "Accepted by every selected step" \
    || fail "--step 1,2 help should show the shared option block" 7
log "Combined help verified for --step 1,2 / 1,3 / 2,3 / 1,2,3."

if [ -f "${ROOT}/mirdeep_p3/cli.py" ]; then
    ( cd "${ROOT}" && python3 -m mirdeep_p3.cli -h >/dev/null 2>&1 ) \
        || fail "python3 -m mirdeep_p3.cli -h failed" 7
    log "Console entry point (mirdeep_p3.cli) verified."
fi

if [ -f "${ROOT}/tests/test_coordinates.py" ]; then
    python3 "${ROOT}/tests/test_coordinates.py" >/dev/null 2>&1 \
        || fail "tests/test_coordinates.py failed" 7
    log "Coordinate unit tests passed (O-05)."
fi
if [ -f "${ROOT}/tests/test_data_index.py" ]; then
    python3 "${ROOT}/tests/test_data_index.py" >/dev/null 2>&1 \
        || fail "tests/test_data_index.py failed" 7
    log "Data-index unit tests passed (O-10)."
fi

if [ -f "${ROOT}/tests/test_mature_boundary.py" ]; then
    python3 "${ROOT}/tests/test_mature_boundary.py" >/dev/null 2>&1 \
        || fail "tests/test_mature_boundary.py failed" 7
    log "Mature-boundary unit tests passed."
fi

INDEX_STATE="$(python3 -c "import sys; sys.path.insert(0, '${ROOT}/src'); from utils.data_index import describe_index, index_status; print(describe_index(index_status('${ROOT}')))" 2>&1)" \
    || fail "data index status could not be read" 7
log "data index: ${INDEX_STATE}"
# The README "Current versions" block and docs/DATA-INDEX.md must quote the live
# versions; a software bump or a new data index that skips them makes the docs lie.
SRC_VERSION="$(sed -n 's/^__version__ *= *"\(.*\)"/\1/p' "${ROOT}/src/version.py" | head -1)"
README_VERSION="$(sed -n 's/^> \*\*Current versions\*\* .*software \*\*v\([^*]*\)\*\*.*/\1/p' "${ROOT}/README.md" | head -1)"
if [ "$README_VERSION" != "$SRC_VERSION" ]; then
    fail "README 'Current versions' says v${README_VERSION:-<none>} but src/version.py says ${SRC_VERSION}; update README.md and docs/DATA-INDEX.md" 7
fi
EXPECTED_TAG="$(sed -n 's/^DATA_INDEX_RELEASE=//p' "${ROOT}/data-index.env" | head -1 | tr -d '[:space:]')"
README_INDEX_TAG="$(sed -n 's/^> \*\*Current versions\*\* .*data index \*\*`\([^`]*\)`\*\*.*/\1/p' "${ROOT}/README.md" | head -1)"
if [ -n "$EXPECTED_TAG" ] && [ "$README_INDEX_TAG" != "$EXPECTED_TAG" ]; then
    fail "README 'Current versions' says index ${README_INDEX_TAG:-<none>} but data-index.env says ${EXPECTED_TAG}" 7
fi
if [ -f "${ROOT}/docs/DATA-INDEX.md" ] && [ -n "$EXPECTED_TAG" ]; then
    grep -q "$EXPECTED_TAG" "${ROOT}/docs/DATA-INDEX.md" \
        || fail "docs/DATA-INDEX.md does not mention the current index tag ${EXPECTED_TAG}" 7
fi
log "Version block in README matches src/version.py (${SRC_VERSION}) and data-index.env (${EXPECTED_TAG})."


# ---- 12. Cleanup -------------------------------------------------------------
if [ "$KEEP_OUTPUT" -eq 0 ]; then
    log "Removing test output (use -k to keep)..."
    rm -rf "$OUT_DIR"
fi

# ---- Done -------------------------------------------------------------------
log "=============================================="
log " ALL TESTS PASSED "
log "=============================================="
exit 0
