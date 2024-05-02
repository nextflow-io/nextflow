#!/bin/bash
set -e 

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

$NXF_CMD run test-complexpaths.nf -c awsbatch.config
[[ -d foo ]] || false
[[ -e 'foo/.alpha' ]] || false
[[ -e 'foo/01_A(R1).fastq' ]] || false
[[ -e 'foo/01_A(R2).fastq' ]] || false
[[ -e 'foo/f1.fa' ]] || false
[[ -e 'foo/f2.fa' ]] || false
[[ -e 'foo/f3.fa' ]] || false
[[ -e 'foo/hello.txt' ]] || false
[[ -e 'foo/sample.html' ]] || false
[[ -e 'foo/sample.zip' ]] || false
[[ -e 'foo/sample_(1 2).vcf' ]] || false

rm -rf foo
$NXF_CMD run test-complexpaths.nf -resume -c awsbatch.config
[[ -d foo ]] || false
[[ -e 'foo/.alpha' ]] || false
[[ -e 'foo/01_A(R1).fastq' ]] || false
[[ -e 'foo/01_A(R2).fastq' ]] || false
[[ -e 'foo/f1.fa' ]] || false
[[ -e 'foo/f2.fa' ]] || false
[[ -e 'foo/f3.fa' ]] || false
[[ -e 'foo/hello.txt' ]] || false
[[ -e 'foo/sample.html' ]] || false
[[ -e 'foo/sample.zip' ]] || false
[[ -e 'foo/sample_(1 2).vcf' ]] || false

$NXF_CMD run test-subdirs.nf -c awsbatch.config

NXF_CLOUDCACHE_PATH=s3://nextflow-ci/cache \
$NXF_CMD run nextflow-io/rnaseq-nf \
    -profile batch \
    -with-report \
    -with-trace \
    -plugins nf-cloudcache
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false

NXF_CLOUDCACHE_PATH=s3://nextflow-ci/cache \
$NXF_CMD run nextflow-io/rnaseq-nf \
    -profile batch \
    -with-report \
    -with-trace \
    -plugins nf-cloudcache \
    -resume
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false
[[ `grep -c 'Cached process > ' .nextflow.log` == 4 ]] || false

## run with fargate + wave
NXF_CLOUDCACHE_PATH=s3://nextflow-ci/cache \
$NXF_CMD run nextflow-io/rnaseq-nf \
    -profile batch \
    -plugins nf-cloudcache,nf-wave \
    -c awsfargate.config

## Test use of job array
NXF_CLOUDCACHE_PATH=s3://nextflow-ci/cache \
$NXF_CMD run nextflow-io/hello \
    -process.array 10 \
    -plugins nf-cloudcache \
    -c awsbatch.config

## Test use of job array using Fusion
$NXF_CMD run nextflow-io/hello \
    -process.array 10 \
    -with-wave \
    -with-fusion \
    -c awsbatch.config