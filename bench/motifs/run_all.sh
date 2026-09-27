#!/bin/bash
# run every motif (or those given), PAR at a time (default 6), each capped at $TMO s
cd "$(dirname "$0")"
OUT=${OUT:-results}; TMO=${TMO:-1200}; PAR=${PAR:-6}; mkdir -p $OUT
MOTIFS=${@:-$(Rscript -e 'source("motifs.R"); cat(names(.motifs))')}
run() { timeout $TMO Rscript run_motif.R $1 $OUT "$EXTRA" > $OUT/$1.log 2>&1 || echo "$1 TIMEOUT/FAIL rc=$?" >> $OUT/$1.log; }
export -f run; export OUT TMO EXTRA
printf "%s\n" $MOTIFS | xargs -P $PAR -I{} bash -c 'run {}'
tail -q -n 1 $OUT/*.log
