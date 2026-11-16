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
# Smoke test of a built nf-agent-pi runner image, needing only Docker: no credentials, no network
# beyond pulling the images. It runs the harness executable the way agent-rpc does, over its JSONL
# protocol, against the scripted model in mock_model.py, which has it `read` a PNG and run `bash`.
# Unit tests run harness/runner.mjs against a stub SDK, so this is what covers the artifact that
# actually ships: the compiled binary and the files it needs at runtime, such as the photon wasm.
#
# Usage: smoke-test.sh <image>
set -euo pipefail

# checked before anything starts: bash 3.2 (macOS /bin/bash) only fails to parse `coproc` when it
# reaches it, after the mock is already running
(( BASH_VERSINFO[0] >= 4 )) || { echo "smoke-test.sh needs bash >= 4 (coproc), found $BASH_VERSION" >&2; exit 1; }

image=${1:?usage: smoke-test.sh <image>}
HERE=$(cd "$(dirname "$0")" && pwd)
MOCK_IMAGE=python:3.12-slim@sha256:e695ab189fa2d082b1c61c69a6b888a5adeca61ceeccad107b1ed2a125721a9e
# must match PiAgentRunner.getLaunchSpec()
HARNESS=/opt/nf-agent-pi/runner
TIMEOUT=120

id=nf-agent-pi-smoke-$$
cleanup() {
    docker rm -f "$id-harness" "$id-mock" >/dev/null 2>&1 || true
    docker network rm "$id" >/dev/null 2>&1 || true
    docker volume rm "$id" >/dev/null 2>&1 || true
}
trap cleanup EXIT

fail() {
    echo "Smoke test failed: $image" >&2
    printf '%s\n' "$@" >&2
    echo "--- mock model logs" >&2
    docker logs "$id-mock" >&2 2>&1 || true
    exit 1
}

docker network create "$id" >/dev/null
docker volume create "$id" >/dev/null
docker run -d --name "$id-mock" --network "$id" -v "$id:/work" -v "$HERE/mock_model.py:/mock_model.py:ro" \
    "$MOCK_IMAGE" python -u /mock_model.py /work >/dev/null
ready=
for _ in $(seq 30); do
    docker exec "$id-mock" python -c 'import socket; socket.create_connection(("localhost", 8080), 1)' 2>/dev/null && { ready=1; break; }
    sleep 1
done
[[ $ready ]] || fail "the mock model did not start accepting connections within 30s"

coproc HARNESS_PROC { docker run -i --rm --name "$id-harness" --network "$id" -v "$id:/work:ro" -w /work "$image" "$HARNESS"; }
# bash unsets the coproc variables once it exits, so keep our own copies
harness_pid=$HARNESS_PROC_PID
exec {harness_out}<&"${HARNESS_PROC[0]}" {harness_in}>&"${HARNESS_PROC[1]}"

echo '{"type":"start","protocolVersion":2,"invocationId":"smoke","apiKey":"smoke-key","spec":{"model":"anthropic/claude-sonnet-4-5","baseUrl":"http://'"$id-mock"':8080","prompt":"smoke test","workDir":"/work","nativeToolNames":["read","bash"],"maxIterations":5}}' >&"$harness_in"

result=
ended=
deadline=$((SECONDS + TIMEOUT))
while [[ -z $result && -z $ended ]]; do
    if (( SECONDS >= deadline )); then
        ended=timeout
    elif IFS= read -r -t $((deadline - SECONDS)) line <&"$harness_out"; then
        case $line in
            *'"type":"complete"'*|*'"type":"error"'*) result=$line ;;
        esac
    elif (( $? > 128 )); then
        ended=timeout
    else
        ended=exited
    fi
done

case $ended in
    timeout) fail "the harness produced no complete or error frame within ${TIMEOUT}s" ;;
    exited)
        status=0
        wait "$harness_pid" || status=$?
        fail "the harness exited (status $status) without a complete or error frame" ;;
esac

if [[ $result == *'"type":"complete"'* && $result == *'"output":"SMOKE-OK"'* ]]; then
    echo "Smoke test passed: $image"
else
    fail "$result"
fi
