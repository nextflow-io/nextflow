set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}


[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
set -x

## pool creation: runs against a pool that doesn't exist yet, in the background and in its
## own launch dir, so the pool's ~6 min start-up overlaps the tests on the shared pool.
## It checks Azure SDK error mapping, which doesn't depend on the JDK, so one job is enough.
if [[ ${TEST_JDK:-17} == 17 ]]; then
  export POOL_ID="nextflow-ci-$(date +%s)-$RANDOM"
  mkdir -p pool-create
  ( cd pool-create && exec timeout -k 5m 20m $NXF_CMD -C ../pool-create.config run ../test-subdirs.nf ) > pool-create.out 2>&1 &
  POOL_CREATE_PID=$!
  # a pool that can't get a node waits forever, hence the timeout; on timeout, or if a foreground
  # test fails, TERM stops the run through Nextflow's shutdown, which deletes the pool
  trap 'kill -TERM $POOL_CREATE_PID 2>/dev/null && { wait $POOL_CREATE_PID; cat pool-create.out; }' EXIT
fi

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

if [[ ${POOL_CREATE_PID:-} ]]; then
  trap - EXIT
  wait $POOL_CREATE_PID || { cat pool-create.out; echo "Pool creation run failed"; exit 1; }
  grep -q "pool with id: $POOL_ID" pool-create/.nextflow.log || { cat pool-create.out; echo "Pool $POOL_ID was not created"; exit 1; }
fi
