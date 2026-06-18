#!/usr/bin/env bash
#
# Build check hook for Nextflow development.
#
# On Stop / SubagentStop, run `make compile` to catch syntax and compilation
# errors. On failure, block and feed a parsed error summary back to Claude.
#
# Reads the Claude Code hook payload as JSON on stdin (requires `jq`).

set -uo pipefail

command -v jq >/dev/null 2>&1 || exit 0

input=$(cat)
event=$(jq -r '.hook_event_name // ""' <<<"$input")

case "$event" in
  Stop|SubagentStop) ;;
  *) exit 0 ;;
esac

start=$SECONDS
if out=$(make compile 2>&1); then
  jq -n --arg s "$(( SECONDS - start ))" \
    '{suppressOutput: true, systemMessage: ("✓ Build check passed (" + $s + "s)")}'
else
  summary=$(printf '%s\n' "$out" | grep -iE 'error|failed|exception' | tail -n 10)
  [[ -n "$summary" ]] || summary=$(printf '%s\n' "$out" | tail -n 20)
  jq -n --arg r "Build check failed:
$summary" \
    '{decision: "block", reason: $r, stopReason: "Build compilation failed"}'
fi
