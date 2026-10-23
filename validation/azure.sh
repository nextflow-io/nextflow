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
## first remove pools and jobs left by runs interrupted before deleting their pool
if command -v az >/dev/null && [[ -n "${AZURE_BATCH_ACCOUNT_KEY:+set}" ]]; then (
  { set +x; } 2>/dev/null
  export AZURE_BATCH_ACCOUNT="$AZURE_BATCH_ACCOUNT_NAME"
  export AZURE_BATCH_ACCESS_KEY="$AZURE_BATCH_ACCOUNT_KEY"
  export AZURE_BATCH_ENDPOINT="https://$AZURE_BATCH_ACCOUNT_NAME.westeurope.batch.azure.com"
  CUTOFF=$(( $(date +%s) - 3*3600 ))
  POOLS=$(timeout 120 az batch pool list --query "[?starts_with(id,'nextflow-ci-')].id" -o tsv) || { echo "WARN: unable to list pools for cleanup"; exit 0; }
  for P in $POOLS; do
    TS=${P#nextflow-ci-}; TS=${TS%%-*}
    [[ $TS =~ ^[0-9]+$ ]] && (( 10#$TS < CUTOFF )) || continue
    for J in $(timeout 120 az batch job list --filter "executionInfo/poolId eq '$P'" --query "[].id" -o tsv); do
      timeout 120 az batch job delete --job-id "$J" --yes || true
    done
    timeout 120 az batch pool delete --pool-id "$P" --yes && echo "Removed stale pool $P" || true
  done
) fi

POOL_ID="nextflow-ci-$(date +%s)-$RANDOM"
POOL_CONFIG=$(mktemp)
trap 'rm -f "$POOL_CONFIG"' EXIT
cat > "$POOL_CONFIG" <<EOF
includeConfig '$PWD/azure.config'
process.queue = '$POOL_ID'
azure {
  batch {
    deletePoolsOnCompletion = true
    pools {
      '$POOL_ID' {
          // a VM family whose core quota the shared nextflow-ci pool doesn't use
          vmType = 'Standard_D2d_v4'
          autoScale = true
      }
    }
  }
}
EOF
$NXF_CMD -C "$POOL_CONFIG" run ./test-subdirs.nf
grep -q "pool with id: $POOL_ID" .nextflow.log || { echo "Pool $POOL_ID was not created"; exit 1; }

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
