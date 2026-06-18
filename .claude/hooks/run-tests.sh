#!/usr/bin/env bash
#
# Test runner hook for Nextflow development.
#
# After an Edit/Write/MultiEdit on a .groovy/.java file under modules/ or
# plugins/, map it to its test class and run it via gradle. A source file runs
# its corresponding *Test class; a test file runs itself. Block on failure.
#
# Reads the Claude Code hook payload as JSON on stdin (requires `jq`).

set -uo pipefail

command -v jq >/dev/null 2>&1 || exit 0

input=$(cat)
tool_name=$(jq -r '.tool_name // ""' <<<"$input")
file_path=$(jq -r '.tool_input.file_path // ""' <<<"$input")

case "$tool_name" in
  Edit|Write|MultiEdit) ;;
  *) exit 0 ;;
esac

[[ -n "$file_path" ]] || exit 0
case "$file_path" in *.groovy|*.java) ;; *) exit 0 ;; esac
case "$file_path" in *modules/*|*plugins/*) ;; *) exit 0 ;; esac

# Resolve the gradle module path.
module=""
case "$file_path" in
  *plugins/nf-*)
    plugin=$(sed -E 's#.*plugins/(nf-[^/]+)/.*#\1#' <<<"$file_path")
    module="plugins:${plugin}" ;;
  *modules/nextflow/*)   module="nextflow" ;;
  *modules/nf-commons/*) module="nf-commons" ;;
  *modules/nf-httpfs/*)  module="nf-httpfs" ;;
  *modules/nf-lang/*)    module="nf-lang" ;;
  *modules/nf-lineage/*) module="nf-lineage" ;;
esac
[[ -n "$module" ]] || exit 0

class=$(basename "$file_path"); class="${class%.*}"
case "$class" in
  *Test) test_class="$class" ;;
  *)     test_class="${class}Test" ;;
esac

if out=$(./gradlew ":${module}:test" --tests "*${test_class}" 2>&1); then
  jq -n --arg f "$(basename "$file_path")" \
    '{suppressOutput: true, systemMessage: ("✓ Tests passed for " + $f)}'
else
  summary=$(printf '%s\n' "$out" | grep -iE 'failed|error|exception|assertion' | tail -n 5)
  [[ -n "$summary" ]] || summary=$(printf '%s\n' "$out" | tail -n 15)
  jq -n --arg f "$(basename "$file_path")" --arg s "$summary" \
    '{decision: "block", reason: ("Tests failed for " + $f + ":\n" + $s)}'
fi
