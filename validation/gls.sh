set -e

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

#
# setup credentials:
#  1. decrypt credentials file
#  2. export required env var
# Note: file was encrypted with the command:
#  gpg --symmetric --cipher-algo AES256 my_secret.json
#
# More details https://help.github.com/en/actions/automating-your-workflow-with-github-actions/creating-and-using-encrypted-secrets
#   
gpg --quiet --batch --yes --decrypt --passphrase=$GOOGLE_SECRET --output google_credentials.json ./google_credentials.gpg
export GOOGLE_APPLICATION_CREDENTIALS=$PWD/google_credentials.json

[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
set -x
$NXF_CMD -C ./gls.config \
    run nextflow-io/rnaseq-nf \
    -with-report \
    -with-trace $OPTS