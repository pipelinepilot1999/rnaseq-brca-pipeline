#!/bin/bash

SAMPLES="SRR1552445 SRR1552447 SRR1552448"
OUTDIR="data/raw"
TMPDIR="data/raw/tmp"

for sample in $SAMPLES; do
    echo "Downloading $sample..."
    fasterq-dump $sample --outdir $OUTDIR --temp $TMPDIR
    echo "Compressing $sample..."
    gzip $OUTDIR/$sample.fastq
    echo "Done: $sample"
done
