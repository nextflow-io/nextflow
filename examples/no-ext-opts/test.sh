#!/bin/bash
# Runs the pipeline twice and asserts the resolved tool args
set -euo pipefail
cd "$(dirname "$0")"

check() { grep -qF -- "$2" "$1" || { echo "FAIL: $1 missing '$2': $(cat "$1")"; exit 1; }; }
reject() { ! grep -qF -- "$2" "$1" || { echo "FAIL: $1 has '$2'"; exit 1; }; }

nextflow run . --input data/samplesheet.csv --fasta data/genome.fa -output-dir out-default -ansi-log false > /dev/null
check out-default/trimgalore/s1.log 'trim_galore --fastqc --clip_r1 10 --paired'   # pipeline default
check out-default/trimgalore/s2.log 'trim_galore --fastqc --clip_r1 5  s2'        # clip_r1 column, SE drops R2 flags
check out-default/trimgalore/s3.log 'trim_galore --custom-trim --paired'          # trimgalore_args column
check out-default/align/s1.align.log 'align --bowtie2 --genome'
test -e out-default/align/s1.sorted.bam                                           # per-call-site prefix
test -e out-default/dedup/s1.dedup.sorted.bam

sleep 5   # global resume: avoid session lock race
nextflow run . --input data/samplesheet.csv --fasta data/genome.fa -output-dir out-override -ansi-log false \
    --rrbs --maxins 500 --clip_r2 3 --skip_dedup --multiqc_title 'My run' --opts.samtools_sort.m=2G > /dev/null
check out-override/trimgalore/s1.log '--rrbs --clip_r1 10 --clip_r2 3'
reject out-override/trimgalore/s2.log '--clip_r2'
check out-override/trimgalore/s3.log '--custom-trim --paired'                     # column beats --args
check out-override/align/s1.align.log '--bowtie2 --maxins 500'
reject out-override/align/s2.align.log '--maxins'
check out-override/align/s1.sorted.bam 'samtools sort -m 2G'
test ! -e out-override/dedup
check out-override/multiqc/multiqc_report.html 'multiqc --title My run'   # stub echo strips the quotes

sleep 5
# per-option merge into defaults: change one, drop one (false), add an unlisted one
nextflow run . --input data/samplesheet.csv --fasta data/genome.fa -output-dir out-opts -ansi-log false \
    --maxins 500 --opts.trimgalore.fastqc=false --opts.trimgalore.quality=30 --opts.align.maxins=700 --opts.align.local=true > /dev/null
check out-opts/trimgalore/s1.log 'trim_galore --clip_r1 10 --quality 30 --paired'
check out-opts/trimgalore/s2.log 'trim_galore --clip_r1 5 --quality 30  s2'
check out-opts/trimgalore/s3.log 'trim_galore --custom-trim --paired'             # column still replaces
check out-opts/align/s1.align.log 'align --bowtie2 --maxins 700 --local --genome'
check out-opts/align/s2.align.log 'align --bowtie2 --maxins 700 --local --genome' # user opts bypass single_end gating

echo PASS
