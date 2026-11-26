#!/usr/bin/env bash
# OLS benchmark: design.matrix's qr() against eight other least squares
# implementations, on four families of designs (numeric only, numeric with
# a factor, one numeric and a factor, numeric x factor interaction) at
# several widths, plus one rank-deficient design. See cell.R for what is
# measured.
#
# Fair: one thread everywhere; design.matrix's na.omit-copy shortcut off;
# one fresh R process per method and design; every result checked against
# lm(). Processes run one at a time on one fast core (this machine is
# hybrid: cores 0, 2, 4, 6 are fast). Several workers on separate fast cores
# were tried and rejected: timings rose 12% to 46% over serial runs, unevenly
# across methods (shared memory bandwidth and lower boost clocks), which
# would bias the comparison. Memory was unaffected.
# Safe: a memory cap and a time limit per process.
#
# Usage: bench/ols/run.sh [cores]   default "0"; one worker per listed core
# Writes bench/ols/results.csv only when every measurement succeeded.
set -u
here=$(cd "$(dirname "$0")" && pwd)
cores=${1:-"0"}
out=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
n=200000
iters=11
designs=()
for w in 5 20 50 100; do designs+=("numeric $w"); done
for w in 25 50 100; do designs+=("mixed $w"); done
for w in 5 20 50 100; do designs+=("factor $w"); done
for w in 6 20 50 100; do designs+=("interaction $w"); done
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
  write.csv(d, args[2], row.names = FALSE)' "$out" "$here/results.csv" || exit 1
Rscript -e 'si <- sessionInfo(); writeLines(c(R.version.string,
  paste("BLAS:", basename(si$BLAS), "(limited to 1 thread)"),
  paste("CPU:", system("lscpu | grep \"Model name\" | sed \"s/.*: *//\"", intern = TRUE)),
  paste("Cores used:", commandArgs(TRUE)[2])), commandArgs(TRUE)[1])' "$here/machine.txt" "$cores"
rm -rf "$out"
echo "done: $got measurements"
