##
##  Copyright 2013-2024, Seqera Labs
##
##  Licensed under the Apache License, Version 2.0 (the "License");
##  you may not use this file except in compliance with the License.
##  You may obtain a copy of the License at
##
##      http://www.apache.org/licenses/LICENSE-2.0
##
##  Unless required by applicable law or agreed to in writing, software
##  distributed under the License is distributed on an "AS IS" BASIS,
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##  See the License for the specific language governing permissions and
##  limitations under the License.
##
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

commands="!{executions}"
IFS=';' read -r -a publications <<< "$commands"
for publication in "${publications[@]}"; do
        eval $publication &
done
wait