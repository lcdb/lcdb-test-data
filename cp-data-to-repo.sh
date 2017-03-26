#!/bin/bash

if [[ -z $1 ]]; then
    echo "usage: $0 BUILT_DATA_PATH"
    exit 1
fi
SRC=$1

# Builds the include file for rsync
cat > .inc <<EOF
+ rnaseq_samples
+ rnaseq_samples/*
+ rnaseq_samples/*/*small*.sorted.bam
+ rnaseq_samples/*/*tiny*.sorted.bam
+ rnaseq_samples/*/*small*fastq.gz
+ rnaseq_samples/*/*tiny*fastq.gz
+ chipseq_samples
+ chipseq_samples/*
+ chipseq_samples/*/*small*.sorted.bam
+ chipseq_samples/*/*tiny*.sorted.bam
+ chipseq_samples/*/*small*fastq.gz
+ chipseq_samples/*/*tiny*fastq.gz
+ seq
+ seq/*small*fa
+ annotation
+ annotation/*small*
- *
EOF

rsync -arvv $SRC --include-from .inc data && rm .inc
