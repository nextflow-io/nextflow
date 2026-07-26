set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

echo "Test smolvm running a task in a microVM"
$NXF_CMD run smolvm.nf -c smolvm.config
