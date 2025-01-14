
set +e
export LC_NUMERIC=C

nxf_publish() {
    local max_retries=$1
    local delay_ms=$2
    local jitter_factor=$3
    local max_delay_ms=$4
    local id=$5
    shift 5
    local attempt=1

    while true; do
        "$@" && echo $id:0 && return # Run the command; if it succeeds, exit the function

        if (( attempt >= max_retries )); then
            echo "Command failed after $attempt attempts." >&2
            echo $id:1
            return
        fi

        # Calculate jitter delay
        local jittered_delay_ms=$(awk -v delay="$delay_ms" -v jitter="$jitter_factor" 'BEGIN { srand(); print ( 1 - jitter + (2 * jitter) * rand() ) * delay }')

        # Ensure the delay does not exceed the maximum delay and convert to seconds for sleep
        local total_delay_sec=$(awk -v delay="$jittered_delay_ms" -v max="$max_delay_ms" 'BEGIN { if (delay > max) delay = max; print delay / 1000 }')
        echo "Command failed. Retrying in $total_delay_sec seconds..." >&2
        sleep "$total_delay_sec"

        # Increment the delay with exponential backoff, clamped to max_delay_ms
        delay_ms=$((delay_ms * 2))
        if (( delay_ms > max_delay_ms )); then
            delay_ms=$max_delay_ms
        fi

        ((attempt++))
    done
}

check_max_parallel_copies() {
    while ((${#pid[@]}>=${max_parallel_copies})); do
        sleep 0.2 2>/dev/null
        local copy=()
        for x in "${pid[@]}"; do
            # if the process exist, keep in the 'copy' array, otherwise wait on it to capture the exit code
            [[ -e /proc/$x ]] && copy+=($x) || wait $x
        done
        pid=("${copy[@]}")
    done
}

commands="!{executions}"
max_parallel_copies="!{max_parallel}"
pid=()
IFS=';' read -r -a publications <<< "$commands"
for publication in "${publications[@]}"; do
    check_max_parallel_copies
    eval $publication &
    pid+=($!)
done
wait