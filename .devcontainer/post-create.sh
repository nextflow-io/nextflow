#!/usr/bin/env bash
set -euo pipefail

echo "[devcontainer] Bootstrapping Nextflow development dependencies..."

# Trigger Gradle wrapper download and prime common local build tasks.
./gradlew --no-daemon --version
./gradlew --no-daemon buildInfo compile

echo "[devcontainer] Done."
