#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/.."

# read the nextflow version
read -r NF_VERSION<VERSION

echo "
-----------------------------------------------
-- Publishing nextflow launcher script to S3 --
-----------------------------------------------
"

# determine publish location based on release type
S3_RELEASE_BUCKET=${S3_RELEASE_BUCKET:-'www2.nextflow.io'}
if [[ "$NF_VERSION" =~ .+(-edge|-EDGE) ]]; then
  S3_RELEASE_DIR="releases/edge"
else
  S3_RELEASE_DIR="releases/latest"
fi

# collect the files we need to publish
files+=(
  'nextflow'
  'nextflow.sha1'
  'nextflow.sha256'
  'nextflow.md5'
  'VERSION'
)

# publish them to S3
for file in "${files[@]}"; do
  filename=$(echo "$file" | tr '[:upper:]' '[:lower:]')
  aws s3 cp "$file" "s3://$S3_RELEASE_BUCKET/$S3_RELEASE_DIR/$filename" \
    --no-progress \
    --storage-class STANDARD \
    --region eu-west-1 \
   --acl public-read
done

echo "Done"
