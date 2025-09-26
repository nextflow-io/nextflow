#!/bin/bash
#
# Nextflow Release Script
#
# This script performs the complete Nextflow release process including:
# - Building and assembling artifacts
# - Uploading to S3 and Maven repositories 
# - Releasing plugins to registry
# - Creating GitHub releases with signed artifacts
#
# REQUIRED SECRETS/ENVIRONMENT VARIABLES FOR GITHUB ACTIONS:
#
# AWS S3 Deployment:
#   NXF_AWS_ACCESS       - AWS Access Key for deploying to s3://www2.nextflow.io
#   NXF_AWS_SECRET       - AWS Secret Key for S3 deployment
#
# Maven Repository (Seqera S3-based):
#   AWS_ACCESS_KEY_ID    - AWS credentials for Maven repository access
#   AWS_SECRET_ACCESS_KEY - AWS secret for Maven repository access
#
# GitHub Integration:
#   GITHUB_TOKEN         - For creating releases and uploading assets
#
# Plugin Registry:
#   NPR_API_URL          - Nextflow Plugin Registry API URL
#   NPR_API_KEY          - Nextflow Plugin Registry API key
#
# Container Registry Authentication:
#   DOCKERHUB_USERNAME   - Docker Hub username for container publishing
#   DOCKERHUB_TOKEN      - Docker Hub token/password for container publishing
#   SEQERA_PUBLIC_CR_PASSWORD - Seqera public container registry password
#
# Usage: Only run when commit message contains '[release]'
#
set -e

echo "=== Starting Nextflow Release Process ==="
echo "Commit message: ${GITHUB_HEAD_COMMIT_MESSAGE:-$(git log -1 --pretty=format:'%s')}"

# Check required environment variables
echo "=== Checking required environment variables ==="
REQUIRED_VARS=(
    "NXF_AWS_ACCESS"
    "NXF_AWS_SECRET"
    "AWS_ACCESS_KEY_ID"
    "AWS_SECRET_ACCESS_KEY"
    "GITHUB_TOKEN"
    "NPR_API_URL"
    "NPR_API_KEY"
)

MISSING_VARS=()
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        MISSING_VARS+=("$var")
    else
        echo "✓ $var is set"
    fi
done

if [ ${#MISSING_VARS[@]} -ne 0 ]; then
    echo "❌ ERROR: The following required environment variables are not set:"
    for var in "${MISSING_VARS[@]}"; do
        echo "  - $var"
    done
    echo "Please ensure all required environment variables are configured before running the release."
    exit 1
fi

echo "✅ All required environment variables are set"

echo "🔧 === Step 1: Assemble, upload, and deploy ==="
make assemble upload deploy
echo "✅ Step 1 completed successfully"
echo ""

echo "🔌 === Step 2: Release plugins ==="
make release-plugins
echo "✅ Step 2 completed successfully"
echo ""

echo "🚀 === Step 3: Final release ==="
make release
echo "✅ Step 3 completed successfully"
echo ""

echo "🎉 === Release process completed successfully ==="
