set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

echo "Test Wave accessing private container repository"
(cd wave-tests/example1; bash run.sh)

echo "Test Wave building a container"
(cd wave-tests/example2; bash run.sh)

echo "Test Wave building from Conda package"
(cd wave-tests/example3; bash run.sh)

echo "Test Wave building from Conda lock file"
(cd wave-tests/example4; bash run.sh)

echo "Test Wave running rnaseq-nf with Fusion on local executor"
(cd wave-tests/example6; bash run.sh)

echo "Test Wave running rnaseq-nf with Fusion on AWS Batch"
(cd wave-tests/example6; bash run-aws.sh)

echo "Test Wave running rnaseq-nf with Fusion on Google Batch"
(cd wave-tests/example6; bash run-gcp.sh)
