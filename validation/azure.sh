set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}


[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
set -x
$NXF_CMD -C ./azure.config \
    run ./test-readspair.nf \
    -with-report \
    -with-trace $OPTS

## complex paths test
$NXF_CMD -C ./azure.config run ./test-complexpaths.nf
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
$NXF_CMD -C ./azure.config run ./test-complexpaths.nf -resume
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

## run test-subdirs inputs/outputs
$NXF_CMD -C ./azure.config run ./test-subdirs.nf

## pool creation: the pool must not exist before the run
POOL_ID="nextflow-ci-$(date +%s)-$RANDOM"
cat > pool-create.config <<EOF
includeConfig 'azure.config'
process.queue = '$POOL_ID'
azure {
  batch {
    deletePoolsOnCompletion = true
    pools {
      '$POOL_ID' {
          vmCount = 1
      }
    }
  }
}
EOF
$NXF_CMD -C ./pool-create.config run ./test-subdirs.nf
[[ `grep -c "pool with id: $POOL_ID" .nextflow.log` == 1 ]] || false

NXF_CLOUDCACHE_PATH=az://my-data/cache \
$NXF_CMD -c ./azure.config \
    run nextflow-io/rnaseq-nf \
    -with-report \
    -with-trace $OPTS \
    -plugins nf-cloudcache
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false

NXF_CLOUDCACHE_PATH=az://my-data/cache \
$NXF_CMD -c ./azure.config \
    run nextflow-io/rnaseq-nf \
    -with-report \
    -with-trace $OPTS \
    -plugins nf-cloudcache \
    -resume
[[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false
[[ `grep -c 'Cached process > ' .nextflow.log` == 4 ]] || false
