set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''

# Each test runs in its own launch dir under $RUNS_DIR, so these paths must be absolute
BASE_DIR=$PWD
CONFIG=$BASE_DIR/google.config
RUNS_DIR=$BASE_DIR/google-runs

test_arrays() {
  $NXF_CMD -C $CONFIG -q run $BASE_DIR/test-arrays.nf > array_output
  [[ `grep 'Hi from the nf-test-array bucket!' -c array_output` == 3 ]] && echo OK || { echo 'Failed array tasks' && false; }
}

test_readspair() {
  $NXF_CMD -C $CONFIG \
      run $BASE_DIR/test-readspair.nf \
      -with-report \
      -with-trace $OPTS
}

check_complexpaths() {
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
}

test_complexpaths() {
  $NXF_CMD -C $CONFIG run $BASE_DIR/test-complexpaths.nf
  check_complexpaths

  rm -rf foo
  $NXF_CMD -C $CONFIG run $BASE_DIR/test-complexpaths.nf -resume
  check_complexpaths
}

test_rnaseq() {
  NXF_CLOUDCACHE_PATH=gs://rnaseq-nf/cache \
  $NXF_CMD -trace nextflow,io.seqera -C $CONFIG \
      run nextflow-io/rnaseq-nf \
      -with-report \
      -with-trace $OPTS \
      -plugins nf-cloudcache
  [[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false

  NXF_CLOUDCACHE_PATH=gs://rnaseq-nf/cache \
  $NXF_CMD -trace nextflow,io.seqera -C $CONFIG \
      run nextflow-io/rnaseq-nf \
      -with-report \
      -with-trace $OPTS \
      -plugins nf-cloudcache \
      -resume
  [[ `grep -c 'Using Nextflow cache factory: nextflow.cache.CloudCacheFactory' .nextflow.log` == 1 ]] || false
  [[ `grep -c 'Cached process > ' .nextflow.log` == 4 ]] || false
}

## Test job array with Fusion
test_hello_fusion() {
  $NXF_CMD -C $CONFIG \
      run nextflow-io/hello \
      -process.array 10 \
      -with-wave \
      -with-fusion
}

## Test job array
test_hello_array() {
  $NXF_CMD -C $CONFIG \
      run nextflow-io/hello \
      -process.array 10
}

TESTS=(test_arrays test_readspair test_complexpaths test_rnaseq test_hello_fusion test_hello_array)

# Pull remote pipelines up front so the concurrent runs don't race on the asset cache
$NXF_CMD pull nextflow-io/rnaseq-nf
$NXF_CMD pull nextflow-io/hello

# Each Batch job spends most of its time waiting for a VM, so run independent tests concurrently
rm -rf $RUNS_DIR
declare -A PIDS
for t in "${TESTS[@]}"; do
  mkdir -p $RUNS_DIR/$t
  ( set -ex; cd $RUNS_DIR/$t; $t ) > $RUNS_DIR/$t/test.log 2>&1 &
  PIDS[$t]=$!
done

FAILED=()
for t in "${TESTS[@]}"; do
  if wait ${PIDS[$t]}; then status=OK; else status=FAILED; FAILED+=("$t"); fi
  echo "==== BEGIN $t ($status) ===="
  cat $RUNS_DIR/$t/test.log
  echo "==== END $t ($status) ===="
done

if (( ${#FAILED[@]} )); then
  for t in "${FAILED[@]}"; do echo "::error::Google Batch test failed: $t (see google-runs/$t/test.log)"; done
  exit 1
fi
