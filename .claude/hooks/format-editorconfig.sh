#!/usr/bin/env bash
#
# EditorConfig enforcement hook for Nextflow development.
#
# After an Edit/Write/MultiEdit, apply the project .editorconfig rules to the
# edited source file with `eclint fix`. Non-blocking: warns on failure, never
# stops Claude. Skips silently if eclint is not installed.
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

[[ -n "$file_path" && -f "$file_path" ]] || exit 0

# Skip dot-directories (.git, .gradle, ...) and build outputs.
case "$file_path" in
  */build/*|*/.*/*) exit 0 ;;
esac

# Only format known source extensions.
case "$file_path" in
  *.groovy|*.java|*.gradle|*.md|*.txt|*.yml|*.yaml|*.json) ;;
  *) exit 0 ;;
esac

# eclint is an optional dependency: skip quietly if it isn't installed.
command -v eclint >/dev/null 2>&1 || exit 0

if out=$(eclint fix "$file_path" 2>&1); then
  jq -n --arg f "$(basename "$file_path")" \
    '{suppressOutput: true, systemMessage: ("✓ EditorConfig formatting applied to " + $f)}'
else
  echo "Warning: eclint failed for $file_path: $out" >&2
  exit 1
fi
