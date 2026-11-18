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
# Dataflow benchmarks. Runs each case with the development launcher, one at a time.
#
#   run.sh [-o OUTDIR] [CASE...]
#
# A CASE is `script:n:heap:passes`, e.g. `dataflow:1000000:256m:time,heap`.
#   script   dataflow | processes | slow-consumer (a .nf file in this directory)
#   n        number of items
#   heap     `default` or a -Xmx value such as `256m`
#   passes   `time` runs without sampling; `heap` also samples the live heap every SAMPLE seconds
#            with `jcmd GC.class_histogram`, which forces a full GC, so its wall time is perturbed
# Without CASEs the default matrix below runs.
#
# Environment: LAUNCH (default: <repo>/launch.sh), TIMEOUT seconds per run (default 600),
# SAMPLE seconds between heap samples (default 2), EXTRA_OPTS appended to NXF_OPTS.
# Every run gets a fresh NXF_HOME, cache and work dir under a `mktemp -d` scratch dir.
# Build first with `make compile`.

set -uo pipefail

BENCH_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_DIR=$(cd "$BENCH_DIR/../../../.." && pwd)
LAUNCH=${LAUNCH:-$REPO_DIR/launch.sh}
TIMEOUT=${TIMEOUT:-600}
SAMPLE=${SAMPLE:-2}
EXTRA_OPTS=${EXTRA_OPTS:-}
OUTDIR=

DEFAULT_CASES=(
    dataflow:1000000:default:time,heap
    dataflow:1000000:256m:time,heap
    processes:10000:default:time,heap
    processes:10000:256m:time,heap
    slow-consumer:1000000:256m:heap
    slow-consumer:1000000:128m:heap
    slow-consumer:30000:default:heap
)

while getopts "o:" opt; do
    case $opt in
        o) OUTDIR=$OPTARG ;;
        *) sed -n '/^# Dataflow benchmarks/,/^# Build first/p' "$0"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
CASES=("$@")
[[ ${#CASES[@]} -eq 0 ]] && CASES=("${DEFAULT_CASES[@]}")
[[ -z $OUTDIR ]] && OUTDIR=$(mktemp -d -t nf-bench-XXXXXX)
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

machine_info() {
    echo "date: $(date -Is)"
    echo "host: $(uname -srm)"
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')"
    echo "nproc: $(nproc)"
    echo "mem: $(free -m | awk '/^Mem:/ {print $2 " MiB total, " $7 " MiB available"}')"
    echo "load: $(cut -d' ' -f1-3 /proc/loadavg)"
    echo "commit: $(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null) $(git -C "$REPO_DIR" status --porcelain 2>/dev/null | wc -l) changed files"
    local java=java
    [[ -n ${JAVA_HOME:-} ]] && java=$JAVA_HOME/bin/java
    "$java" -version 2>&1 | sed 's/^/java: /'
}

# seconds from the `Elapsed (wall clock)` line of /usr/bin/time -v
wall_seconds() {
    awk -F': ' '/Elapsed \(wall clock\)/ {
        n = split($2, p, ":"); s = 0
        for( i = 1; i <= n; i++ ) s = s * 60 + p[i]
        printf "%.1f", s }' "$1"
}

# count, full count, total pause ms, and max heap after GC (MiB) from a -Xlog:gc log
gc_summary() {
    awk '/Pause/ {
            c++; if( $0 ~ /Pause Full/ ) f++
            ms = $NF; sub(/ms$/, "", ms); t += ms
            if( match($0, /->[0-9]+M/) ) { a = substr($0, RSTART + 2, RLENGTH - 3) + 0; if( a > m ) m = a }
        }
        END { printf "%d\t%d\t%.0f\t%d", c, f, t, m }' "$1"
}

# sample the live heap of the JVM tagged with $1 until it exits; writes "elapsed_s live_bytes" lines
sample_heap() {
    local tag=$1 out=$2 pid= t0=$SECONDS max=0
    for _ in $(seq 1 200); do
        pid=$(pgrep -n -f -- "-Dbench.tag=$tag") && break
        sleep 0.05
    done
    [[ -z $pid ]] && return
    while kill -0 "$pid" 2>/dev/null; do
        sleep "$SAMPLE"
        local hist total
        hist=$(jcmd "$pid" GC.class_histogram 2>/dev/null) || continue
        total=$(awk '/^Total/ {print $3}' <<< "$hist")
        [[ -z $total ]] && continue
        echo "$((SECONDS - t0)) $total" >> "$out/heap-samples.txt"
        if (( total > max )); then
            max=$total
            head -n 23 <<< "$hist" > "$out/heap-peak-histogram.txt"
        fi
    done
}

run_case() {
    local spec=$1 script n heap passes
    IFS=: read -r script n heap passes <<< "$spec"
    for pass in ${passes//,/ }; do
        local label="${script}-n${n}-${heap}-${pass}"
        local out="$OUTDIR/$label"
        local scratch tag
        mkdir -p "$out"
        scratch=$(mktemp -d -t nf-bench-run-XXXXXX)
        tag=$(basename "$scratch")
        local opts="-Dbench.tag=$tag -XX:+ExitOnOutOfMemoryError -Xlog:gc:file=$out/gc.log"
        [[ $heap != default ]] && opts="-Xmx$heap $opts"
        [[ -n $EXTRA_OPTS ]] && opts="$opts $EXTRA_OPTS"
        local cmd=("$LAUNCH" run "$BENCH_DIR/$script.nf" --n "$n")
        echo "cd $scratch && NXF_HOME=$scratch/home NXF_CACHE_DIR=$scratch/cache NXF_WORK=$scratch/work NXF_OPTS='$opts' ${cmd[*]}" > "$out/command.txt"
        echo ">> $label" >&2

        (
            cd "$scratch" || exit 1
            export NXF_HOME=$scratch/home NXF_CACHE_DIR=$scratch/cache NXF_WORK=$scratch/work NXF_OPTS=$opts NXF_AGENT_MODE=false
            /usr/bin/time -v -o "$out/time.txt" timeout -k 30 "$TIMEOUT" "${cmd[@]}" > "$out/stdout.txt" 2>&1
        ) &
        local runner=$!
        [[ $pass == heap ]] && sample_heap "$tag" "$out"
        wait $runner

        local status wall rss gc live result oom=no
        status=$(awk -F': ' '/Exit status/ {print $2}' "$out/time.txt")
        wall=$(wall_seconds "$out/time.txt")
        rss=$(awk -F': ' '/Maximum resident set size/ {printf "%d", $2 / 1024}' "$out/time.txt")
        gc=$(gc_summary "$out/gc.log")
        live=-
        [[ -s $out/heap-samples.txt ]] && live=$(sort -k2 -n "$out/heap-samples.txt" | tail -1 | awk '{printf "%d", $2 / 1048576}')
        grep -q OutOfMemoryError "$out/stdout.txt" && oom=yes
        [[ $status == 124 || $status == 137 ]] && oom=timeout
        result=$(grep -o 'result=.*' "$out/stdout.txt" | head -1)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$script" "$n" "$heap" "$pass" "$status" "$oom" "$wall" "$rss" "$live" "$gc" "${result:--}" "$out" >> "$OUTDIR/results.tsv"
        tail -1 "$OUTDIR/results.tsv" >&2
    done
}

machine_info > "$OUTDIR/machine.txt"
printf 'script\tn\theap\tpass\texit\toom\twall_s\tpeak_rss_mib\tlive_heap_mib\tgc_count\tfull_gc\tgc_pause_ms\tmax_heap_after_gc_mib\tresult\tdir\n' > "$OUTDIR/results.tsv"
for c in "${CASES[@]}"; do
    run_case "$c"
done
echo "results: $OUTDIR/results.tsv" >&2
