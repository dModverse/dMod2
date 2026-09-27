#!/bin/bash
# run every motif (or those given) in parallel, each capped at $TMO seconds
cd "$(dirname "$0")"
OUT=${OUT:-results}; TMO=${TMO:-1200}; mkdir -p $OUT
MOTIFS=${@:-$(Rscript -e 'source("motifs.R"); cat(names(.motifs))')}
for m in $MOTIFS; do
  ( timeout $TMO Rscript run_motif.R $m $OUT "$EXTRA" > $OUT/$m.log 2>&1 || echo "$m TIMEOUT/FAIL rc=$?" >> $OUT/$m.log ) &
done
wait
tail -q -n 1 $OUT/*.log
