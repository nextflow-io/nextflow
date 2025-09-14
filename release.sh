#!/bin/bash
set -e
set -x

echo "=== Starting Nextflow Release Process ==="
echo "Commit message: ${GITHUB_HEAD_COMMIT_MESSAGE:-$(git log -1 --pretty=format:'%s')}"

# Verify we're in release mode
if [[ ! "${GITHUB_HEAD_COMMIT_MESSAGE:-$(git log -1 --pretty=format:'%s')}" =~ \[release\] ]]; then
    echo "ERROR: This script should only be run when commit message contains '[release]'"
    exit 1
fi

echo "=== Step 1: Clean, assemble, upload, and deploy ==="
make clean assemble upload deploy

echo "=== Step 2: Release plugins ==="
make release-plugins

echo "=== Step 3: Final release ==="
make release

echo "=== Release process completed successfully ==="