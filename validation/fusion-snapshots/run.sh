#!/bin/bash
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
# Runs the snapshot test, checks the counter survived the reclaim, then checks -resume caches it.
set -euo pipefail
cd "$(dirname "$0")"
NXF=${NXF:-nextflow}
STEPS=${STEPS:-600}

for v in E2E_FIS_TEMPLATE E2E_WORK_DIR E2E_SPOT_QUEUE E2E_JOB_ROLE AWS_REGION TOWER_ACCESS_TOKEN; do
    [ -n "${!v:-}" ] || { echo "missing env: $v"; exit 1; }
done

fail() { echo "FAIL: $*"; exit 1; }

$NXF run main.nf -ansi-log false --steps $STEPS

# a restart from scratch writes 1 again; a restore may repeat the few lines written after the snapshot
[ "$(grep -c '^1 ' results/count.txt)" -eq 1 ] || fail "counter started twice, task restarted instead of restoring"
[ "$(tail -1 results/count.txt | cut -d' ' -f1)" -eq $STEPS ] || fail "counter did not reach $STEPS"
[ "$(cut -d' ' -f2 results/count.txt | uniq | wc -l)" -ge 2 ] \
    || fail "task never moved to a new instance, the interruption did not happen"
echo "steps per instance (a restore continues the count on the next instance):"
awk '$2!=prev{if(prev)print "  " first "-" last " " prev; first=$1; prev=$2} {last=$1} END{print "  " first "-" last " " prev}' results/count.txt
grep -q 'COUNT.*COMPLETED.*1$' trace.txt || fail "expected one Nextflow attempt: $(cat trace.txt)"

$NXF run main.nf -ansi-log false --steps $STEPS -resume
grep -q 'COUNT.*CACHED' trace.txt || fail "COUNT not cached on resume: $(cat trace.txt)"

echo "PASS"
