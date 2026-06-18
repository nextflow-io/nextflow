#!/usr/bin/env bash
#
# IntelliJ IDEA formatter hook for Nextflow development.
#
# After an Edit/Write/MultiEdit on a .groovy file, run IntelliJ IDEA's headless
# formatter to match project code style. IDEA is located via the IDEA_SH
# environment variable or common install paths.
#
# Non-blocking: skips silently if IDEA is not installed (optional), and only
# warns (never blocks) if the formatter errors. Wired with `async: true`, so it
# runs in the background and does not hold up the conversation.
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

case "$file_path" in *.groovy) ;; *) exit 0 ;; esac
[[ -f "$file_path" ]] || exit 0
case "$file_path" in */build/*|*/.*/*) exit 0 ;; esac

find_idea() {
  if [[ -n "${IDEA_SH:-}" && -x "$IDEA_SH" ]]; then
    echo "$IDEA_SH"; return 0
  fi
  local candidates=(
    "/Applications/IntelliJ IDEA.app/Contents/MacOS/idea"
    "/Applications/IntelliJ IDEA CE.app/Contents/MacOS/idea"
    "/Applications/IntelliJ IDEA Ultimate.app/Contents/MacOS/idea"
    "/usr/local/bin/idea.sh"
    "/opt/idea/bin/idea.sh"
    "/snap/intellij-idea-community/current/bin/idea.sh"
    "/snap/intellij-idea-ultimate/current/bin/idea.sh"
    "$HOME/.local/share/JetBrains/Toolbox/apps/IDEA-U/bin/idea.sh"
    "$HOME/.local/share/JetBrains/Toolbox/apps/IDEA-C/bin/idea.sh"
  )
  local c
  for c in "${candidates[@]}"; do
    [[ -x "$c" ]] && { echo "$c"; return 0; }
  done
  command -v idea 2>/dev/null && return 0
  return 1
}

# IDEA is an optional dependency: skip quietly if it isn't installed.
idea=$(find_idea) || exit 0

# Prefer a sibling format.sh, else fall back to `idea format`.
idea_dir=$(dirname "$idea")
if [[ -x "$idea_dir/format.sh" ]]; then
  fmt=("$idea_dir/format.sh" -allowDefaults "$file_path")
else
  fmt=("$idea" format -allowDefaults "$file_path")
fi

if ! out=$("${fmt[@]}" 2>&1); then
  # Non-blocking: a formatting failure should never stop the workflow.
  echo "Warning: IDEA formatter failed for $file_path: ${out:-(no output)}" >&2
  exit 1
fi
