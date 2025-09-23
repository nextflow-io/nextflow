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

## Verify we're in release mode
#if [[ ! "${GITHUB_HEAD_COMMIT_MESSAGE:-$(git log -1 --pretty=format:'%s')}" =~ \[release\] ]]; then
#    echo "ERROR: This script should only be run when commit message contains '[release]'"
#    exit 1
#fi

echo "=== Step 1: Assemble, upload, and deploy ==="
make assemble upload deploy

echo "=== Step 2: Release plugins ==="
make release-plugins

echo "=== Step 3: Final release ==="
make release

echo "=== Release process completed successfully ==="
