#!/bin/bash
set -e 

S3_REPO_BASE=$1

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}


export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

## Test S3 git remote integration: push and run
echo "Testing S3 git remote integration"

# Create timestamp for unique test repo
TIMESTAMP=$(date +%s)
TEST_REPO="test-repo-${TIMESTAMP}"
S3_REPO="$S3_REPO_BASE/${TEST_REPO}"
TEMP_DIR=$(mktemp -d -t nf-s3-test-XXXXXX)
## Test S3 git remote integration: push and run
echo "Testing S3 git remote integration in ${TEMP_DIR}"

remove() {
    echo "Removing ${S3_REPO}"
    $NXF_CMD fs rm "${S3_REPO}"
    echo "Removing ${TEMP_DIR}"
    rm -rf "${TEMP_DIR}"
}

# Copy test pipeline to temp directory
cp -r s3-remote-test-repo "${TEMP_DIR}/"

# Copy config
cp awsbatch.config "${TEMP_DIR}/nextflow.config"

# define trap to remove when exit
trap remove EXIT

cd ${TEMP_DIR}/s3-remote-test-repo
echo "Pushing pipeline to ${S3_REPO} (with explicit repo URL)"
$NXF_CMD push "${S3_REPO}" -r main -c -y

echo "Running pipeline from S3 remote"
$NXF_CMD -q run "${S3_REPO}" -r main | tee stdout1
[[ `grep -c "Hey! Bonjour world!" stdout1` == 1 ]] || false
[[ `grep -c "Hey! Ciao world!" stdout1` == 1 ]] || false
[[ `grep -c "Hey! Hello world!" stdout1` == 1 ]] || false
[[ `grep -c "Hey! Hola world!" stdout1` == 1 ]] || false

echo "Modifying pipeline message"
sed -i "s/Hey!/Hey there!/g" ${TEMP_DIR}/s3-remote-test-repo/main.nf

echo "Pushing modified pipeline to ${S3_REPO} (auto-detect from git remote)"
$NXF_CMD push -m "Update greeting message" -c -y

echo "Running modified pipeline from S3 remote"
$NXF_CMD -q run "${S3_REPO}" -r main -latest | tee stdout2
[[ `grep -c "Hey there! Bonjour world!" stdout2` == 1 ]] || false
[[ `grep -c "Hey there! Ciao world!" stdout2` == 1 ]] || false
[[ `grep -c "Hey there! Hello world!" stdout2` == 1 ]] || false
[[ `grep -c "Hey there! Hola world!" stdout2` == 1 ]] || false

echo "S3 git remote integration test completed successfully"


