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
{ set +x; } 2>/dev/null
if command -v az >/dev/null && [[ -n "${AZURE_BATCH_ACCOUNT_KEY:-}" ]]; then
  BATCH_OPTS=(--account-name "$AZURE_BATCH_ACCOUNT_NAME" --account-key "$AZURE_BATCH_ACCOUNT_KEY" --account-endpoint "https://$AZURE_BATCH_ACCOUNT_NAME.westeurope.batch.azure.com")
  CUTOFF=$(( $(date +%s) - 3*3600 ))
  for P in $(az batch pool list "${BATCH_OPTS[@]}" --query "[?starts_with(id,'nextflow-ci-')].id" -o tsv 2>/dev/null); do
    TS=${P#nextflow-ci-}; TS=${TS%%-*}
    if [[ $TS =~ ^[0-9]+$ ]] && (( TS < CUTOFF )); then
      for J in $(az batch job list "${BATCH_OPTS[@]}" --query "[?poolInfo.poolId=='$P'].id" -o tsv 2>/dev/null); do
        az batch job delete "${BATCH_OPTS[@]}" --job-id "$J" --yes || true
      done
      az batch pool delete "${BATCH_OPTS[@]}" --pool-id "$P" --yes || true
      echo "Removed stale pool $P"
    fi
  done
fi
set -x

POOL_ID="nextflow-ci-$(date +%s)-$RANDOM"
POOL_CONFIG=$(mktemp)
cat > $POOL_CONFIG <<EOF
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
$NXF_CMD -C $POOL_CONFIG run ./test-subdirs.nf
rm -f $POOL_CONFIG
grep -q "pool with id: $POOL_ID" .nextflow.log

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
