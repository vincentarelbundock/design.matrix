#!/usr/bin/env bash
# OLS benchmark: design.matrix's qr() against other least squares
# implementations, from formula and data and from a prebuilt matrix.
#
# Fair: one thread everywhere; design.matrix's na.omit-copy shortcut off;
# one fresh R process per method and design (memory = peak resident memory
# of the first call, all allocations included; then the repeats are timed);
# every result checked against lm().
# Safe: a memory cap and a time limit per process; designs whose dense
# model matrix would exceed 2 GB are skipped.
#
# Usage: bench/ols/run.sh   (from the package root, with the package and
# the competitors installed). Writes bench/ols/results.csv, only when every
# expected measurement succeeded; otherwise exits 1 and leaves it as is.
set -u
here=$(cd "$(dirname "$0")" && pwd)
out=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
n=200000
declare -A P=([p3]=3 [p5]=5 [p21]=21 [p40]=40 [p119]=119 [p401]=401)
declare -A NM=([formula]=9 [prebuilt]=8)
failed=0
expected=0
for des in p3 p5 p21 p40 p119 p401; do
  bytes=$(( n * ${P[$des]} * 8 ))
  if [ $bytes -gt 2000000000 ]; then echo "skip $des"; continue; fi
  expected=$(( expected + NM[formula] + NM[prebuilt] ))
  ( ulimit -v 25000000; timeout 1200 Rscript "$here/ref.R" $des $n "$out" ) || { echo "fail ref $des"; failed=1; continue; }
  it=5; [ $des = p401 ] && it=3
  for mode in formula prebuilt; do
    for mi in $(seq 1 ${NM[$mode]}); do
      ( ulimit -v 25000000; timeout 1200 Rscript "$here/cell.R" $des $n $it $mode "$out" $mi ) \
        > /dev/null 2>&1 || { echo "fail $mode $des method $mi"; failed=1; }
    done
  done
done
got=$(ls "$out" | grep -c -E '^(formula|prebuilt).*csv$')
if [ $failed -ne 0 ] || [ "$got" -ne "$expected" ]; then
  echo "incomplete: $got of $expected measurements; results.csv not updated (raw files in $out)" >&2
  exit 1
fi
# A method that stops with an error is a result (recorded with its message);
# a row without timing and without a message is a harness failure.
Rscript -e 'args <- commandArgs(TRUE); f <- list.files(args[1], "^(formula|prebuilt).*csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(f, read.csv))
  bad <- is.na(d$ms) & (is.na(d$note) | d$note == "")
  if (any(bad)) { print(d[bad, ]); quit(status = 1) }
  err <- is.na(d$ms)
  if (any(err)) { message("methods that stopped with an error:"); print(d[err, c("mode", "design", "method", "note")]) }
  d <- d[order(d$mode, d$p, d$ms), ]
  write.csv(d, args[2], row.names = FALSE)' "$out" "$here/results.csv" || exit 1
Rscript -e 'si <- sessionInfo(); writeLines(c(R.version.string, paste("BLAS:", basename(si$BLAS)), paste("CPU:", system("lscpu | grep \"Model name\" | sed \"s/.*: *//\"", intern = TRUE))), commandArgs(TRUE)[1])' "$here/machine.txt"
rm -rf "$out"
