#!/usr/bin/env bash
#
# Copyright 2013-2026, Seqera Labs
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Requires Docker, kind, kubectl and jq. Run make compile in the checkout first.
# Usage: bash validation/k8s-automatic-retries.sh [launcher] [separate|baseline]
# A temporary kind cluster is deleted on exit. The logs are kept for inspection.

set -euo pipefail

fail() { printf '%s\n' "$*" >&2; exit 1; }
root=$(cd "$(dirname "$0")/.." && pwd -P)
launcher=${1:-$root/launch.sh}
expected=${2:-separate}
[[ $# -le 2 && ( $expected == separate || $expected == baseline ) ]] || fail "Expected separate or baseline"
[[ -x $launcher ]] || fail "Launcher is not executable: $launcher"
launcher="$(cd "$(dirname "$launcher")" && pwd -P)/$(basename "$launcher")"
for tool in docker kind kubectl jq; do
    command -v "$tool" >/dev/null || fail "Missing command: $tool"
done

run_dir=$(mktemp -d "${TMPDIR:-/tmp}/nf-kind-retries.XXXXXX")
run_dir=$(cd "$run_dir" && pwd -P)
cluster=$(basename "$run_dir" | tr '[:upper:].' '[:lower:]-')
export KUBECONFIG="$run_dir/kubeconfig"
workflow_pid=
cleanup() {
    if [[ -n $workflow_pid ]] && kill -0 "$workflow_pid" 2>/dev/null; then
        kill "$workflow_pid" 2>/dev/null || true
        wait "$workflow_pid" 2>/dev/null || true
    fi
    kind delete cluster --name "$cluster" --kubeconfig "$KUBECONFIG"
    rm -f "$KUBECONFIG"
    printf 'Logs and task output: %s\n' "$run_dir"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
printf 'Creating %s for the %s check\n' "$cluster" "$expected"

shared="$run_dir/shared"
mkdir -p "$shared"
cat > "$run_dir/kind.yaml" <<EOF
kind: Cluster
apiVersion: kind.x-k8s.io/v1alpha4
nodes:
- role: control-plane
  extraMounts:
  - hostPath: "$shared"
    containerPath: "$shared"
EOF
kind create cluster --name "$cluster" --kubeconfig "$KUBECONFIG" --config "$run_dir/kind.yaml" --wait 120s

context="kind-$cluster"
namespace=retry-test
kube() { kubectl --context "$context" --namespace "$namespace" --request-timeout=15s "$@"; }
kube create namespace "$namespace"
kube apply -f - <<EOF
apiVersion: v1
kind: PersistentVolume
metadata:
  name: retry-test
spec:
  capacity:
    storage: 1Gi
  accessModes: [ReadWriteOnce]
  storageClassName: ""
  hostPath:
    path: "$shared"
---
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: retry-test
spec:
  accessModes: [ReadWriteOnce]
  storageClassName: ""
  volumeName: retry-test
  resources:
    requests:
      storage: 1Gi
EOF

cat > "$run_dir/nextflow.config" <<EOF
params.shared = '$shared'
process {
    executor = 'k8s'
    container = 'ubuntu:24.04'
    errorStrategy = 'retry'
    maxRetries = 3
}
k8s {
    context = '$context'
    namespace = '$namespace'
    storageClaimName = 'retry-test'
    storageMountPath = '$shared'
    computeResourceType = 'Pod'
}
EOF
cat > "$run_dir/main.nf" <<'EOF'
process CHECK_RETRY {
    script:
    """
    if [ '${task.attempt}' -eq 1 ]; then
        while [ ! -f '${params.shared}/release-'\$(hostname) ]; do sleep 1; done
    fi
    echo 'application-attempt=${task.attempt}' >> '${params.shared}/attempts'
    if [ '${task.attempt}' -eq 1 ]; then exit 42; fi
    """
}

workflow {
    CHECK_RETRY()
}
EOF

(
    cd "$run_dir"
    export NXF_ANSI_LOG=false
    exec "$launcher" -log nextflow.log -C nextflow.config run main.nf -work-dir "$shared/work"
) > "$run_dir/stdout.log" 2>&1 &
workflow_pid=$!

wait_for_pod() {
    local deadline=$((SECONDS + 180)) pod
    while :; do
        kill -0 "$workflow_pid" 2>/dev/null || fail "Nextflow exited before starting the next pod"
        (( SECONDS < deadline )) || fail "Timed out waiting for the next pod"
        pod=$(kube get pods -l nextflow.io/app=nextflow --field-selector=status.phase=Running -o json |
            jq -r 'first(.items[] | select(.metadata.deletionTimestamp == null) | .metadata.name) // empty')
        if [[ -n $pod ]]; then
            printf '%s\n' "$pod"
            return
        fi
        sleep 1
    done
}

# A zero-grace eviction can remove the pod from the API before its container stops.
# Release only the fourth pod so an evicted container cannot run the command.
for attempt in 1 2 3; do
    pod=$(wait_for_pod)
    printf 'Evicting pod %s of 3 (%s)\n' "$attempt" "$pod"
    jq -n --arg pod "$pod" --arg namespace "$namespace" '{
        apiVersion: "policy/v1",
        kind: "Eviction",
        metadata: {name: $pod, namespace: $namespace},
        deleteOptions: {gracePeriodSeconds: 0}
    }' | kube create --raw "/api/v1/namespaces/$namespace/pods/$pod/eviction" -f - >/dev/null
    deadline=$((SECONDS + 30))
    while :; do
        remaining=$(kube get pod "$pod" --ignore-not-found -o name)
        [[ -z $remaining ]] && break
        (( SECONDS < deadline )) || fail "Timed out waiting for pod $pod to disappear"
        sleep 1
    done
done
pod=$(wait_for_pod)
touch "$shared/release-$pod"

deadline=$((SECONDS + 180))
while kill -0 "$workflow_pid" 2>/dev/null; do
    (( SECONDS < deadline )) || fail "Timed out waiting for the workflow to finish"
    sleep 1
done
status=0
wait "$workflow_pid" || status=$?
workflow_pid=

log="$run_dir/nextflow.log"
retries=$(grep -c 'Execution is retried$' "$log" || true)
ordinary=$(grep -Ec 'Execution is retried \([0-9]+\)$' "$log" || true)
[[ $retries == 3 ]] || fail "Expected 3 automatic retries, found $retries"
[[ -f $shared/attempts ]] || fail "The application command did not run after the evictions"
attempts=$(cat "$shared/attempts")
case $expected in
    separate)
        [[ $status == 0 && $ordinary == 1 && $attempts == $'application-attempt=1\napplication-attempt=2' ]] ||
            fail "Expected one application retry at attempt 2 and a successful run (status=$status, ordinary retries=$ordinary)"
        ;;
    baseline)
        task_exit=$(awk '/^Command exit status:$/ { getline; print $1; exit }' "$run_dir/stdout.log")
        [[ $status == 1 && $task_exit == 42 && $ordinary == 0 && $attempts == application-attempt=1 ]] ||
            fail "Expected the application to fail at attempt 1 with exit 42 (status=$status, task exit=$task_exit, ordinary retries=$ordinary)"
        ;;
esac
printf 'Verified %s: %s automatic retries, %s ordinary retries, Nextflow exit %s\n' \
    "$expected" "$retries" "$ordinary" "$status"
printf '%s\n' "$attempts"
