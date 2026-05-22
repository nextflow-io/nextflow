#!/usr/bin/env bash
set -euo pipefail

# Contract: this file is valid Bash as written. Nextflow must not rewrite
# ${...} expressions in the source before execution.
#
# Proposed runtime API:
#   nextflow_input    associative array for process inputs
#   nextflow_output   associative array for process outputs
#   nextflow_params   associative array for process params
#   nextflow          associative array for task metadata/resources

input_fa="${nextflow_input[reads]}"
output_fa="${nextflow_output[sorted]}"
context_txt="${nextflow_output[context]}"

awk '/^>/{print; next} {print toupper($0)}' "$input_fa" > "$output_fa"

cat > "$context_txt" <<EOF
language=bash
reads=$input_fa
sorted=$output_fa
cpus=${nextflow[cpus]}
attempt=${nextflow[attempt]}
EOF
