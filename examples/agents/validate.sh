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
# End-to-end validation of examples/agents/* against the local development
# build, on the local executor and/or the Kubernetes profile.
#
#   make compile assemble                     # .launch.classpath + plugin dev libs
#   ./examples/agents/validate.sh             # every example, both modes
#   ./examples/agents/validate.sh -m local 01_structured-output 04_tool
#   ./examples/agents/validate.sh -r          # ... and check -resume replays from cache
#
# Requires OPENAI_API_KEY (the examples call a real model) and, in local mode, a
# running Docker daemon: the `pi` runner ships no host-local runtime, so every
# agent task runs in the runner image the nf-agent-pi plugin declares -- nothing
# needs Node on this machine. The k8s
# mode additionally needs the cluster, Wave and S3 access that
# examples/agents/k8s-local.config describes.
#
# `-r` adds the companion check: a second run, with -resume, in the SAME directory,
# which must replay from cache rather than call the model again. It costs one extra
# run per example rather than two, because the fresh run above IS the first half of
# the pair.
#
# Exits non-zero if any run fails, so it can gate a release check.
set -u

MODES="local k8s"
TIMEOUT=2400
JOBS=3
DRY=0
RESUME=0
RESULTS=${AGENT_VALIDATION_DIR:-}

usage() {
  cat <<'TXT'
Validate examples/agents/* against the local development build.

  usage: ./examples/agents/validate.sh [-m local|k8s|both] [-t secs] [-j n] [-o dir] [-r] [-n] [example ...]

    -m  mode; default both
    -t  per-run timeout in seconds; default 2400
    -j  runs in parallel; default 3
    -o  results directory; default build/agent-validation
    -r  also re-run each example with -resume and check it replays from cache
    -n  dry run: run every setup check and print the plan, launch nothing

Run 'make compile assemble' first. Requires OPENAI_API_KEY; local mode needs a
running Docker daemon (every agent task runs in the pi runner image the plugin
declares), and k8s mode also needs the cluster, Wave and
S3 access examples/agents/k8s-local.config expects.

A full run calls a real model once per agent invocation and can take an hour, so
-n first: it answers "is this machine set up, and is the example list what I
think it is" for free. -r roughly doubles the wall clock.
TXT
  exit "${1:-0}"
}

while getopts ':m:t:j:o:rnh' opt; do
  case $opt in
    m) case $OPTARG in local|k8s) MODES=$OPTARG ;; both) MODES="local k8s" ;;
         *) echo "Unknown mode: $OPTARG (use local, k8s or both)" >&2; exit 2 ;; esac ;;
    t) TIMEOUT=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    o) RESULTS=$OPTARG ;;
    r) RESUME=1 ;;
    n) DRY=1 ;;
    h) usage 0 ;;
    *) usage 2 ;;
  esac
done
shift $((OPTIND-1))

AGENTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE_DIR=$(cd "$AGENTS_DIR/../.." && pwd)
RESULTS=${RESULTS:-$BASE_DIR/build/agent-validation}

# `timeout` is GNU coreutils; on macOS it arrives as gtimeout with brew.
TIMEOUT_BIN=$(command -v timeout || command -v gtimeout) || {
  echo "ERROR: no 'timeout' command found (macOS: brew install coreutils)" >&2; exit 2; }

[[ -f $BASE_DIR/.launch.classpath ]] || {
  echo "ERROR: missing .launch.classpath -- run 'make compile' first" >&2; exit 2; }
[[ -d $BASE_DIR/plugins/nf-agent-pi/build/target/libs ]] || {
  echo "ERROR: nf-agent-pi is not built -- run 'make assemble' first" >&2; exit 2; }
[[ -n ${OPENAI_API_KEY:-} ]] || {
  echo "ERROR: OPENAI_API_KEY is not set -- the examples call a real model" >&2; exit 2; }
# The `pi` runner carries no host-local runtime, so every agent task runs in the runner
# image through the engine its config enables -- Docker, for these examples. Check it here
# rather than let each run fail on the first task.
[[ $MODES != *local* ]] || docker info >/dev/null 2>&1 || {
  echo "ERROR: the Docker daemon is not available -- local mode runs every agent task" >&2
  echo "       in the pi runner image the nf-agent-pi plugin declares" >&2; exit 2; }

# An example is a directory holding a `main.nf` -- NOT merely a directory. Discovering by
# directory alone also picks up this script's own results dir, a dot-dir left behind by some
# other tool, and a half-built example carrying only `data/`; each of those then fails every
# single run with `Missing project main script`, which reads like a product bug in the report.
EXAMPLES=("$@")
if [[ ${#EXAMPLES[@]} -eq 0 ]]; then
  while IFS= read -r d; do EXAMPLES+=("$(basename "$d")"); done \
    < <(find "$AGENTS_DIR" -mindepth 1 -maxdepth 1 -type d -exec test -f '{}/main.nf' \; -print | sort)
fi

# Discovery can no longer yield one of these, so this only ever catches a name typed on the
# command line. Say so before burning a run on it.
notexamples=()
for ex in "${EXAMPLES[@]}"; do
  [[ -f $AGENTS_DIR/$ex/main.nf ]] || notexamples+=("$ex")
done
if [[ ${#notexamples[@]} -gt 0 ]]; then
  echo "ERROR: not an agent example (no main.nf): ${notexamples[*]}" >&2
  exit 2
fi

# Four examples assemble a real genome and need an input FASTQ that is gitignored on purpose
# (examples/data/README.md). Without it they fail on a missing input long before reaching any
# agent, which reads like a product bug.
#
# Checked through each example's `data` symlink rather than at $DATA_DIR directly: the symlink
# is what the run actually resolves, so a missing or broken one has to fail here too -- and it
# is the only failure mode this layout adds over a per-example copy.
FASTQ_URL=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz
# ${AGENTS_DIR%/*} rather than $AGENTS_DIR/..: AGENTS_DIR is already absolute (line 87), and
# this keeps the paths printed below readable instead of ending in `/agents/../data`.
DATA_DIR=${AGENTS_DIR%/*}/data
missing=()
for ex in "${EXAMPLES[@]}"; do
  case $ex in
    11_contig-filter)  [[ -f $AGENTS_DIR/$ex/data/sample.fastq.gz ]] || missing+=("$ex") ;;
    07_module-as-tool|09_goal-directed|12_isolate-triage)
                       [[ -f $AGENTS_DIR/$ex/data/sample.fastq ]]    || missing+=("$ex") ;;
  esac
done
if [[ ${#missing[@]} -gt 0 ]]; then
  echo "ERROR: missing input FASTQ for: ${missing[*]}" >&2
  echo "One fetch serves all of them -- they share $DATA_DIR through a 'data' symlink." >&2
  echo "Both forms are needed (11_contig-filter reads the gzipped one):" >&2
  echo "  mkdir -p $DATA_DIR" >&2
  echo "  curl -sL $FASTQ_URL -o $DATA_DIR/sample.fastq.gz" >&2
  echo "  gunzip -c $DATA_DIR/sample.fastq.gz > $DATA_DIR/sample.fastq" >&2
  for ex in "${missing[@]}"; do
    [[ -L $AGENTS_DIR/$ex/data ]] \
      || echo "NOTE: $ex/data is not a symlink -- expected one pointing at ../../data" >&2
  done
  exit 2
fi

# One run, in its own directory: concurrent runs must not share .nextflow.log,
# .nextflow/cache or work/. With -r the pair runs SEQUENTIALLY inside this function --
# the resume run must see the cache the fresh run just wrote -- while different
# examples still run concurrently.
run_one() {
  local ex=$1 mode=$2 dir=$RESULTS/$mode/$1 t0=$SECONDS
  rm -rf "$dir" && mkdir -p "$dir" || return 2
  local args=( run -ansi-log false "$AGENTS_DIR/$ex" )
  [[ $mode == k8s ]] && args+=( -c "$AGENTS_DIR/k8s-local.config" )
  ( cd "$dir" && "$TIMEOUT_BIN" -s TERM -k 30 "$TIMEOUT" "$BASE_DIR/launch.sh" "${args[@]}" ) \
    > "$dir/console.log" 2>&1
  echo $? > "$dir/status"
  echo $((SECONDS-t0)) > "$dir/elapsed"
  # Snapshot the fresh run's log: a resume run would rotate it to .nextflow.log.1 and the
  # transport signals below must keep describing the run that actually called the model.
  cp "$dir/.nextflow.log" "$dir/fresh.nextflow.log" 2>/dev/null
  (( RESUME )) || return 0
  t0=$SECONDS
  ( cd "$dir" && "$TIMEOUT_BIN" -s TERM -k 30 "$TIMEOUT" "$BASE_DIR/launch.sh" "${args[@]}" -resume ) \
    > "$dir/resume.log" 2>&1
  echo $? > "$dir/resume-status"
  echo $((SECONDS-t0)) > "$dir/resume-elapsed"
}

# Every check above has now passed, so a dry run has already answered the question it
# exists for -- is this machine set up, and is the example list the expected one. Print
# the plan and stop before the first model call.
if (( DRY )); then
  for mode in $MODES; do
    echo "== $mode: ${#EXAMPLES[@]} examples, $JOBS at a time (dry run)$( (( RESUME )) && echo ', + -resume pass')"
    for ex in "${EXAMPLES[@]}"; do
      cfg=""; [[ $mode == k8s ]] && cfg=" -c $AGENTS_DIR/k8s-local.config"
      printf '   %-24s %s\n' "$ex" "launch.sh run -ansi-log false $AGENTS_DIR/$ex$cfg"
    done
  done
  echo
  echo "dry run: setup checks passed; nothing launched"
  echo "results would go to: $RESULTS"
  exit 0
fi

for mode in $MODES; do
  echo "== $mode: ${#EXAMPLES[@]} examples, $JOBS at a time$( (( RESUME )) && echo ' (fresh + resume)')"
  for ex in "${EXAMPLES[@]}"; do
    while [[ $(jobs -rp | wc -l) -ge $JOBS ]]; do wait -n 2>/dev/null || sleep 2; done
    run_one "$ex" "$mode" &
  done
  wait
done

# -- report. Every example runs through the plugin-hosted RPC broker, so the
#    transport signals are as much the result as the exit code. `remote` is now
#    always true -- the task is always containerized -- so a `remote=false` means a
#    stale build resolved a host-local launch path that no longer exists.
printf '\n%-24s' EXAMPLE; for m in $MODES; do printf '%-34s' "$(echo "$m" | tr a-z A-Z)"; done; echo
failed=0
for ex in "${EXAMPLES[@]}"; do
  printf '%-24s' "$ex"
  for mode in $MODES; do
    dir=$RESULTS/$mode/$ex log=$RESULTS/$mode/$ex/fresh.nextflow.log
    rc=$(cat "$dir/status" 2>/dev/null || echo '?')
    secs=$(cat "$dir/elapsed" 2>/dev/null || echo '?')
    # `grep -c` prints the count AND exits 1 on no match, so reset rather than
    # appending a fallback -- `|| echo 0` yields a two-line "0\n0" that compares unequal
    reg=$(grep -c 'Registering agent RPC invocation' "$log" 2>/dev/null) || reg=0
    rej=$(grep -c 'Rejected agent RPC connection' "$log" 2>/dev/null) || rej=0
    remote=$(grep -o 'remote=[a-z]*' "$log" 2>/dev/null | sort -u | sed 's/remote=//' | paste -sd, -)
    [[ $rc == 0 && $rej == 0 ]] || failed=$((failed+1))
    printf '%-34s' "$([[ $rc == 0 ]] && echo ok || echo "FAIL($rc)") ${secs}s reg=$reg remote=${remote:--} rej=$rej"
  done
  echo
done

# -- resume report.
#
# The signal that the model was NOT called again is `completed=0`, NOT reg=0. A capability
# is minted while the task script is GENERATED, which happens before `storeDir` and the
# resume cache are consulted, so a cache-hit task still registers one and leaves it
# unconsumed -- the broker says exactly that at shutdown. Gating on reg would fail every
# correct resume.
if (( RESUME )); then
  # What the workflow itself printed, with Nextflow's chrome removed. `view` output is the
  # only thing that must match across the pair; staging/pull lines legitimately appear on the
  # fresh run and not on the resumed one.
  answer_of() {
    perl -pe 's/(\[(?:WARN|ERROR|INFO|PIPELINE|WORKDIR|PROCESS|SUCCESS|FAILED)\b)/\n$1/g' "$1" 2>/dev/null \
      | grep -vE '^\[(WARN|ERROR|INFO|PIPELINE|WORKDIR|PROCESS|SUCCESS|FAILED)\b' \
      | grep -vE '^(N E X T F L O W|Launching|executor >|Staging foreign file:|Pulling |Uploading |Downloading |Monitor)' \
      | sed 's/[[:space:]]*$//' | grep -vE '^$'
  }

  printf '\n%-24s' 'EXAMPLE (-resume)'; for m in $MODES; do printf '%-34s' "$(echo "$m" | tr a-z A-Z)"; done; echo
  for ex in "${EXAMPLES[@]}"; do
    printf '%-24s' "$ex"
    for mode in $MODES; do
      dir=$RESULTS/$mode/$ex
      rc2=$(cat "$dir/resume-status" 2>/dev/null || echo '?')
      sum=$(grep -o 'completed=[0-9]* failed=[0-9]* cached=[0-9]*' "$dir/resume.log" 2>/dev/null | tail -1)
      c2=$(sed -E 's/completed=([0-9]+).*/\1/' <<<"$sum"); k2=$(sed -E 's/.*cached=([0-9]+).*/\1/' <<<"$sum")
      answer_of "$dir/console.log" > "$dir/answer-fresh.txt"
      answer_of "$dir/resume.log"  > "$dir/answer-resume.txt"
      # Independent agents' view() output interleaves in a nondeterministic ORDER, so compare
      # order-insensitively. And when nothing executes there are no progress lines to break
      # two view outputs apart, so one whose text lacks a trailing newline abuts the next --
      # same content, different layout. Label that `same*` rather than hide it.
      if diff -q <(sort "$dir/answer-fresh.txt") <(sort "$dir/answer-resume.txt") >/dev/null 2>&1; then
        same=same
      elif [[ "$(tr -d '[:space:]' < "$dir/answer-fresh.txt")" == "$(tr -d '[:space:]' < "$dir/answer-resume.txt")" ]]; then
        same='same*'
      else
        same=DIFFER
      fi
      ok=1
      [[ $rc2 == 0 ]]                       || ok=0
      [[ -n $sum && ${c2:-1} == 0 ]]        || ok=0   # nothing may re-execute
      [[ ${k2:-0} -gt 0 ]]                  || ok=0   # something must come from cache
      [[ $same == same || $same == 'same*' ]] || ok=0
      (( ok )) || failed=$((failed+1))
      printf '%-34s' "$( (( ok )) && echo ok || echo FAIL) completed=${c2:-?} cached=${k2:-?} $same"
    done
    echo
  done
  echo
  echo "  completed=0 is the no-model-call signal; 'same*' means the answers match ignoring whitespace"
fi

# The broker moved from nextflow core into nf-agent-pi; the old core logger
# reappearing means a stale class is being resolved.
stale=$(grep -rl 'DEBUG nextflow\.agent\.AgentRpcBroker' "$RESULTS" 2>/dev/null | wc -l | tr -d ' ')
[[ $stale == 0 ]] || { echo "WARNING: $stale run(s) logged the pre-move core AgentRpcBroker"; failed=$((failed+1)); }

echo
echo "results: $RESULTS"
if [[ $failed -gt 0 ]]; then echo "FAILED: $failed"; exit 1; fi
echo "all runs passed"
