set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
set -x

$NXF_CMD -C ./google.config -q run ./test-arrays.nf > array_output
[[ `grep 'Hi from the nf-test-array bucket!' -c array_output` == 3 ]] && echo OK || { echo 'Failed array tasks' && false; }


$NXF_CMD -C ./google.config \
    run ./test-readspair.nf \
    -with-report \
    -with-trace $OPTS

## complex paths test
$NXF_CMD -C ./google.config run ./test-complexpaths.nf
# validate
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
$NXF_CMD -C ./google.config run ./test-complexpaths.nf -resume
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

NXF_CLOUDCACHE_PATH=gs://rnaseq-nf/cache \
$NXF_CMD -trace nextflow,io.seqera -C ./google.config \
    run nextflow-io/rnaseq-nf \
    -with-report \
    -with-trace $OPTS \
    -plugins nf-cloudcache
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false

NXF_CLOUDCACHE_PATH=gs://rnaseq-nf/cache \
$NXF_CMD -trace nextflow,io.seqera -C ./google.config \
    run nextflow-io/rnaseq-nf \
    -with-report \
    -with-trace $OPTS \
    -plugins nf-cloudcache \
    -resume
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false
[[ `grep -c 'Cached process > ' .nextflow.log` == 4 ]] || false

## Test job array with Fusion
$NXF_CMD -C ./google.config \
    run nextflow-io/hello \
    -process.array 10 \
    -with-wave \
    -with-fusion
    
## Test job array
$NXF_CMD -C ./google.config \
    run nextflow-io/hello \
    -process.array 10


