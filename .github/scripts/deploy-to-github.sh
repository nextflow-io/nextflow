#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/../.."

# read the nextflow version
read -r NF_VERSION<VERSION

function get_change_notes() {
  # extract the most recent changelog notes from 'changelog.txt',
  # which should have been manually updated for the release
  #
  # TODO eventually replace with proper automated changelog generation
  is_relevant=false
  while IFS="" read -r line || [ -n "$line" ]
  do
    # ignore all lines until we find the correct version
    if [[ "$line" == "$NF_VERSION "* ]]; then
      is_relevant=true
    fi

    # then, if line starts with a dash, add to notes
    if [[ "$is_relevant" == true && "$line" == -* ]] ; then
      printf '%s\n' "$line"
    fi
    # until the first empty line
    if [[ "$is_relevant" == true && -z "$line" ]]; then
      break
    fi
  done < changelog.txt
}

echo "Publishing nextflow release to github"

# create a github (pre)release and attach launcher and dist files
# use --verify-tag to fail if tag doesn't exist
notes=$(get_change_notes)

# if edge version, mark as pre-release
if [[ "$NF_VERSION" =~ .+(-edge|-EDGE) ]]; then
  prerelease='--prerelease'
else
  prerelease=
fi

gh release create \
  --draft $prerelease \
  --title "Version $NF_VERSION" \
  --notes "$notes" \
  --verify-tag \
  "v$NF_VERSION" \
  nextflow \
  "build/releases/nextflow-$NF_VERSION-dist"

echo "Done"
