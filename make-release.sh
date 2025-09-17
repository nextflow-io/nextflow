#!/usr/bin/env bash
set -e

# -----------------------------------------------------------------------------
# Nextflow release entrypoint
#
# This script starts a Nextflow release, and is executed by the `make release`
# command. It will help guide you through the manual steps required to perform
# a release, and then trigger a release build on the CI system (Github Actions).
# -----------------------------------------------------------------------------

cd "$(dirname "$0")"

# read the nextflow version
read -r NF_VERSION<VERSION

# read the plugin versions
plugins=()
plugin_versions=()
function get_plugin_version() {
  grep 'Plugin-Version' "$1/src/resources/META-INF/MANIFEST.MF" |\
    cut -d ':' -f 2 |\
    xargs
}
for plugin in plugins/nf-*; do
  if [[ -d "$plugin" ]]; then
    plugins+=("$plugin")
    plugin_versions+=("$(get_plugin_version "$plugin")")
  fi
done

# prompt user to check versions are correct
echo "
------------------------
This will release:

  nextflow: $NF_VERSION
"
for i in "${!plugins[@]}"; do
  plugin_name=$(basename "${plugins[$i]}")
  echo "  $plugin_name: ${plugin_versions[$i]}"
done
echo "------------------------

If these are not the versions you want to release, type 'no' and update the below
files (without committing), then run this script again:
  .
  ├── VERSION
  ├── nextflow
  ├── changelog.txt
  ├── plugins/nf-*/
  │   ├── src/main/resources/META-INF/MANIFEST.MF
  │   └── changelog.txt
  └── modules/nextflow/
      └── src/main/resources/META-INF/plugins-info.txt
"

echo -n "Type 'yes' to proceed: "
read -r proceed
if [[ "$proceed" != "yes" ]]; then
  echo "Aborting"
  exit
fi

# update the digest files before committing anything
./gradlew makeDigest

# create the release commit, including text '[release]' to trigger the github action
# the github action workflow will perform the following tasks
# - tag the release
# - build the release artifacts
# - deploy the release artifacts to maven/docker/S3/etc
# - create a (pre-release) github release
# - deploy the plugins, and update the plugin index
git commit -s -am "Release $NF_VERSION [release]"
# push the current branch is to the remote origin
git push origin HEAD

echo "
------------------------------------------------------------
Release commit pushed.
This should trigger a github release workflow.
Once the workflow is complete, you should check and publish
the draft github release.
------------------------------------------------------------
"
