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

        # Calculate jitter as a factor of the delay
        local min_factor=$(bc <<< "1 - $jitter_factor")
        local max_factor=$(bc <<< "1 + $jitter_factor")
        local random_factor=$(awk -v min="$min_factor" -v max="$max_factor" 'BEGIN{srand(); print min + (max-min)*rand()}')

        # Apply the jitter factor to the delay
        local jittered_delay_ms=$(bc <<< "$delay_ms * $random_factor")

        # Ensure the delay does not exceed the maximum delay
        if (( $(bc <<< "$jittered_delay_ms > $max_delay_ms") == 1 )); then
            jittered_delay_ms=$max_delay_ms
        fi

        # Convert milliseconds to seconds for sleep
        local total_delay_sec=$(bc <<< "scale=3; $jittered_delay_ms / 1000")

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

IFS=$';'
publications=('!{executions}')
for publication in "${publications[@]}"; do
        eval $publication &
done
unset IFS
wait