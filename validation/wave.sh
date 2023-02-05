set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

(cd wave-tests/example1; bash run.sh)
(cd wave-tests/example2; bash run.sh)
(cd wave-tests/example6; bash run.sh)
