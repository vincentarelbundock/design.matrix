#!/usr/bin/env bash
# OLS benchmark: xtx's qr() against eight other least squares
# implementations, on four families of designs (numeric only, numeric with
# a factor, one numeric and a factor, numeric x factor interaction) at
# several widths, plus one rank-deficient design. See cell.R for what is
# measured.
#
# Fair: one thread everywhere; xtx's na.omit-copy shortcut off;
# one fresh R process per method and design; every result checked against
# lm(). Processes run one at a time on one fast core (this machine is
# hybrid: cores 0, 2, 4, 6 are fast). Several workers on separate fast cores
# were tried and rejected: timings rose 12% to 46% over serial runs, unevenly
# across methods (shared memory bandwidth and lower boost clocks), which
# would bias the comparison. Memory was unaffected.
# Safe: a memory cap and a time limit per process.
#
# Usage: bench/ols/run.sh [cores]   default "0"; one worker per listed core
# Writes bench/ols/results.csv (and machine.txt) only when every measurement
# succeeded.
#
# BENCH_SIZE=quick (make bench-quick): two widths per family and 3 timed
# calls per method, to check the scripts. Quick results never go to
# bench/ols/: they go to $BENCH_OUT/ols (default
# ${TMPDIR:-/tmp}/xtx-bench-quick/ols), printed at the end.
set -u
here=$(cd "$(dirname "$0")" && pwd)
cores=${1:-"0"}
size=${BENCH_SIZE:-full}
out=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Rows: 50,000 by default (OLS_BENCH_N): the slowest implementations take
# seconds per fit at the widest designs, and each fit is timed several times
# in a fresh process, so larger data mostly add waiting.
n=${OLS_BENCH_N:-50000}
case "$size" in
  full)
    dest=$here
    iters=11
    w_numeric="5 20 50 100"; w_mixed="25 50 100"; w_factor="5 20 50 100"; w_interaction="6 20 50 100" ;;
  quick)
    dest=${BENCH_OUT:-${TMPDIR:-/tmp}/xtx-bench-quick}/ols
    mkdir -p "$dest" || exit 1
    iters=3
    w_numeric="5 50"; w_mixed="25 100"; w_factor="5 50"; w_interaction="6 50" ;;
  *)
    echo "BENCH_SIZE must be full or quick, not $size" >&2
    exit 1 ;;
esac
designs=()
for w in $w_numeric; do designs+=("numeric $w"); done
for w in $w_mixed; do designs+=("mixed $w"); done
for w in $w_factor; do designs+=("factor $w"); done
for w in $w_interaction; do designs+=("interaction $w"); done
designs+=("aliased 49")

jobs_file="$out/jobs"
: > "$jobs_file"
for d in "${designs[@]}"; do echo "ref $d" >> "$jobs_file"; done
for d in "${designs[@]}"; do
  for mi in 1 2 3 4 5 6 7 8 9; do echo "cell $d formula $mi" >> "$jobs_file"; done
  for mi in 1 2 3 4 5 6 7 8; do echo "cell $d prebuilt $mi" >> "$jobs_file"; done
done
expected=$(grep -c '^cell' "$jobs_file")

run_job() {  # core, then the job's words
  local core=$1; shift
  if [ "$1" = ref ]; then
    ( ulimit -v 25000000; timeout 1200 taskset -c "$core" Rscript "$here/ref.R" "$2" "$3" $n "$out" ) \
      > /dev/null 2>&1 || echo "fail: $*"
  else
    ( ulimit -v 25000000; timeout 1200 taskset -c "$core" Rscript "$here/cell.R" "$2" "$3" $n $iters "$4" "$out" "$5" ) \
      > /dev/null 2>&1 || echo "fail: $*"
  fi
}

# one queue per core, jobs dealt round-robin; references first
run_phase() {  # pattern
  local i=0 pids=() c
  read -ra cs <<< "$cores"
  for c in "${cs[@]}"; do
    ( k=0; grep "^$1" "$jobs_file" | while read -r job; do
        [ $(( k % ${#cs[@]} )) -eq $i ] && run_job "$c" $job
        k=$(( k + 1 ))
      done ) &
    pids+=($!); i=$(( i + 1 ))
  done
  wait "${pids[@]}"
}
run_phase ref | tee "$out/failures"
run_phase cell | tee -a "$out/failures"

got=$(ls "$out" | grep -c -E '^(formula|prebuilt)_.*csv$')
if [ -s "$out/failures" ] || [ "$got" -ne "$expected" ]; then
  echo "incomplete: $got of $expected measurements; results.csv not updated (raw files in $out)" >&2
  exit 1
fi
# A method that stops with an error is a result (recorded with its message);
# a row without timing and without a message is a harness failure.
Rscript -e 'args <- commandArgs(TRUE); f <- list.files(args[1], "^(formula|prebuilt)_.*csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(f, read.csv))
  bad <- is.na(d$ms) & (is.na(d$note) | d$note == "")
  if (any(bad)) { print(d[bad, ]); quit(status = 1) }
  err <- is.na(d$ms)
  if (any(err)) { message("methods that stopped with an error:"); print(d[err, c("mode", "family", "width", "method", "note")]) }
  d <- d[order(d$mode, d$family, d$width, d$ms), ]
  write.csv(d, args[2], row.names = FALSE)' "$out" "$dest/results.csv" || exit 1
Rscript -e 'si <- sessionInfo(); writeLines(c(R.version.string,
  paste("BLAS:", basename(si$BLAS), "(limited to 1 thread)"),
  paste("CPU:", system("lscpu | grep \"Model name\" | sed \"s/.*: *//\"", intern = TRUE)),
  paste("Cores used:", commandArgs(TRUE)[2]),
  if (commandArgs(TRUE)[3] == "quick") "Quick run (BENCH_SIZE=quick): fewer widths, 3 timed calls"),
  commandArgs(TRUE)[1])' "$dest/machine.txt" "$cores" "$size"
rm -rf "$out"
echo "done: $got measurements in $dest/results.csv"
