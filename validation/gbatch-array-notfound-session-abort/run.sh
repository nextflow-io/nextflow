#!/usr/bin/env bash
# usage: ./run.sh <gcp-project> <gs://work-bucket/path> [location]
set -euo pipefail
PROJECT=${1:?gcp project}
WORK=${2:?gs:// work directory}
LOCATION=${3:-us-central1}
${NXF_CMD:-nextflow} run demo.nf -w "$WORK" --project "$PROJECT" --location "$LOCATION"
