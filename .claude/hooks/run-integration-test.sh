#!/usr/bin/env bash
#
# Integration test hook for Nextflow development.
#
# When an integration test case (`tests/<name>.nf`) is edited, run just that one
# case through the existing integration runner (`tests/checks/run.sh`) instead of
# the whole ~115-case suite.
#
# Wired with `asyncRewake: true`: the run happens in the background (so it never
# blocks the conversation) and, on failure (exit code 2), wakes Claude with the
# failure output. The `if` filter in settings.json already scopes this to
# `tests/*.nf` edits; the checks below repeat that scoping because the `if`
# filter is best-effort and fails open.
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

# Only integration cases: a *.nf file directly under tests/.
[[ -n "$file_path" ]] || exit 0
case "$file_path" in *.nf) ;; *) exit 0 ;; esac
[[ "$(basename "$(dirname "$file_path")")" == "tests" ]] || exit 0

root="${CLAUDE_PROJECT_DIR:-$(git -C "$(dirname "$file_path")" rev-parse --show-toplevel 2>/dev/null)}"
[[ -n "$root" && -f "$root/tests/checks/run.sh" ]] || exit 0

name=$(basename "$file_path")

# Use the development launcher; enable Docker only when it is available
# (run.sh skips Docker-only cases when WITH_DOCKER is empty).
export NXF_CMD="$root/launch.sh"
command -v docker >/dev/null 2>&1 && export WITH_DOCKER="-with-docker"

out=$(cd "$root/tests/checks" && bash run.sh "../$name" 2>&1)
status=$?

if [[ $status -ne 0 ]]; then
  {
    echo "Integration test failed: $name"
    printf '%s\n' "$out" | tail -n 30
  } >&2
  exit 2   # asyncRewake: wake Claude with the stderr above
fi
exit 0
