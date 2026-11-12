#!/bin/bash
#
# Deletes the `nextflow-ci-<epoch>-*` pools, and their jobs, that azure.sh runs
# leave behind when they are killed before Nextflow can delete the pool.
# Leaked pools count against the Batch account's pool quota.
#
# Set DRY_RUN=1 to only list what would be deleted.
#
set -euo pipefail

command -v az >/dev/null || { echo "ERROR: Azure CLI (az) not found"; exit 1; }
: "${AZURE_BATCH_ACCOUNT_NAME:?}"
: "${AZURE_BATCH_ACCOUNT_KEY:?}"

export AZURE_BATCH_ACCOUNT="$AZURE_BATCH_ACCOUNT_NAME"
export AZURE_BATCH_ACCESS_KEY="$AZURE_BATCH_ACCOUNT_KEY"
export AZURE_BATCH_ENDPOINT="https://$AZURE_BATCH_ACCOUNT_NAME.westeurope.batch.azure.com"

# pools younger than this may belong to a run still in progress
CUTOFF=$(( $(date +%s) - 3*3600 ))

remove_pool() {
  local J JOBS
  JOBS=$(az batch job list --filter "executionInfo/poolId eq '$1' and state eq 'active'" --query "[].id" -o tsv) || return 1
  for J in $JOBS; do
    az batch job delete --job-id "$J" --yes || return 1
  done
  az batch pool delete --pool-id "$1" --yes
}

POOLS=$(az batch pool list --query "[?starts_with(id,'nextflow-ci-') && state=='active'].id" -o tsv)
STALE=0
FAILED=0
for P in $POOLS; do
  TS=${P#nextflow-ci-}; TS=${TS%%-*}
  [[ $TS =~ ^[0-9]+$ ]] && (( 10#$TS < CUTOFF )) || continue
  STALE=$((STALE+1))
  if [[ ${DRY_RUN:-} ]]; then
    echo "Would remove stale pool $P"
  elif remove_pool "$P"; then
    echo "Removed stale pool $P"
  else
    echo "ERROR: unable to remove stale pool $P"; FAILED=$((FAILED+1))
  fi
done
echo "Stale pools found: $STALE, failed to remove: $FAILED"
(( FAILED == 0 ))
