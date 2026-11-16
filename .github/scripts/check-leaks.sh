#!/bin/bash
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

#
# Scan CI test reports and artifacts for a leaked canary value.
#
# The CI sets NXF_CANARY_TOKEN to a random per-run value in the environment
# of the test steps. If any file about to be uploaded contains the canary,
# a test has dumped the environment (and thus the real secrets next to it).
#
# Usage: check-leaks.sh PATH...
#
# Each PATH can be a directory, a regular file or a tar archive (*.tar,
# *.tar.gz, *.tgz). Missing paths are skipped with a warning.
#
# The canary is searched raw and base64-encoded (standard and URL-safe
# alphabets, at the three possible byte alignments, so that e.g. a Basic
# auth header `user:<canary>` is also detected). The canary must only use
# [A-Za-z0-9_-] so that URL, JSON and HTML escaping leave it unchanged.
#
# Only file paths and match counts are printed, never the matching content.
#
# Exit status: 0 = no leak, 1 = leak detected, 2 = usage error.
#
set -uo pipefail

CANARY_VAR=NXF_CANARY_TOKEN
CANARY=${!CANARY_VAR:-}

if [ -z "$CANARY" ]; then
  echo "check-leaks: $CANARY_VAR is not set" >&2
  exit 2
fi
if [[ ! "$CANARY" =~ ^[A-Za-z0-9_-]{16,}$ ]]; then
  echo "check-leaks: $CANARY_VAR must be at least 16 chars of [A-Za-z0-9_-]" >&2
  exit 2
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
PATTERNS="$WORK/patterns"

# Base64 of the canary preceded by 0, 1 or 2 arbitrary bytes. Drop the
# leading group (it depends on the preceding bytes) when there is a prefix,
# and the trailing partial group (it depends on the following bytes).
{
  printf '%s\n' "$CANARY"
  for prefix in "" "x" "xx"; do
    encoded=$(printf '%s%s' "$prefix" "$CANARY" | base64 -w0)
    # shellcheck disable=SC2017 # truncation to whole groups is intended
    full=$(( (${#prefix} + ${#CANARY}) / 3 * 4 ))
    skip=$(( ${#prefix} > 0 ? 4 : 0 ))
    core=${encoded:skip:full-skip}
    printf '%s\n' "$core"
    printf '%s\n' "$core" | tr '+/' '-_'
  done
} | sort -u > "$PATTERNS"

LEAKS=0

# scan_dir DIR LABEL: report files under DIR containing any pattern,
# displaying paths relative to LABEL
scan_dir() {
  local dir=$1 label=$2 file count
  while IFS= read -r -d '' file; do
    count=$(grep -caF -f "$PATTERNS" -- "$file")
    echo "LEAK: ${label}${file#"$dir"} ($count matching lines)"
    LEAKS=$((LEAKS + 1))
  done < <(grep -rlaZF -f "$PATTERNS" -- "$dir" 2>/dev/null)
}

for path in "$@"; do
  if [ -d "$path" ]; then
    scan_dir "$path" "$path"
  elif [ -f "$path" ] && [[ "$path" =~ \.(tar|tar\.gz|tgz)$ ]]; then
    dest=$(mktemp -d -p "$WORK")
    # `tar -xf` detects compression, CI archives are not always gzipped
    if ! tar -xf "$path" -C "$dest" --no-same-owner --no-same-permissions 2>/dev/null; then
      echo "WARN: failed to fully extract $path, scanning partial content" >&2
    fi
    chmod -R u+rX "$dest" 2>/dev/null
    scan_dir "$dest" "$path!"
    rm -rf "$dest"
  elif [ -f "$path" ]; then
    if grep -qaF -f "$PATTERNS" -- "$path"; then
      echo "LEAK: $path ($(grep -caF -f "$PATTERNS" -- "$path") matching lines)"
      LEAKS=$((LEAKS + 1))
    fi
  else
    echo "WARN: $path not found, skipping" >&2
  fi
done

if [ "$LEAKS" -gt 0 ]; then
  echo "check-leaks: canary found in $LEAKS file(s)"
  exit 1
fi
echo "check-leaks: no leak found"
