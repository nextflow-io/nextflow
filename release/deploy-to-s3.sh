#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/.."

# read the nextflow version
read -r NF_VERSION<VERSION

echo "
-------------------------------------------
-- Uploading nextflow distribution to S3 --
-------------------------------------------
"

S3_RELEASE_BUCKET=${S3_RELEASE_BUCKET:-'www2.nextflow.io'}
S3_RELEASE_DIR="releases/v$NF_VERSION"

# check if the release already exists
release_exists=false
aws s3api head-object --bucket "$S3_RELEASE_BUCKET" --key "$S3_RELEASE_DIR/nextflow" > /dev/null 2>&1 \
  && release_exists=true

if [[ $release_exists == true ]]; then
  echo "Version $NF_VERSION already deployed to S3, skipping"
  exit
fi

# collect files to deploy
files=(build/releases/nextflow-"$NF_VERSION"-*)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR - can't find any files to upload"
  exit 1
fi
files+=(
  'nextflow'
  'nextflow.sha1'
  'nextflow.sha256'
  'nextflow.md5'
)

# upload them to s3 bucket
for file in "${files[@]}"; do
  filename=$(basename "$file")
  aws s3 cp "$file" "s3://$S3_RELEASE_BUCKET/$S3_RELEASE_DIR/$filename" \
    --no-progress \
    --storage-class STANDARD \
    --region eu-west-1 \
    --acl public-read
done

echo "Done"
