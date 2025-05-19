#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/.."

GH_ORG=${GH_ORG:-'nextflow-io'}

# function to extract plugin version from manifest
function get_plugin_version() {
  grep 'Plugin-Version' "$1/src/resources/META-INF/MANIFEST.MF" |\
    cut -d ':' -f 2 |\
    xargs
}

# deploy plugin artifacts to github releases
echo "
----------------------------------
-- Publishing plugins to github --
----------------------------------
"

for plugin in plugins/nf-*; do
  if [[ -d "$plugin" ]]; then
    # get plugin name and version
    plugin_name=$(basename "$plugin")
    plugin_repo="$GH_ORG/$plugin_name"
    plugin_version=$(get_plugin_version "$plugin")

    # check if release already exists
    release_exists=false
    gh release view --repo "$plugin_repo" "$plugin_version" > /dev/null 2>&1 \
      && release_exists=true

    # if not exists, create github release, with zip & meta json files
    if [[ $release_exists == true ]]; then
      echo "Plugin $plugin_name $plugin_version already deployed to github, skipping"
    else
      gh release create \
        --repo "$plugin_repo" \
        --title "Version $plugin_version" \
        "$plugin_version" \
        "$plugin/build/libs/$plugin_name-$plugin_version.zip" \
        "$plugin/build/libs/$plugin_name-$plugin_version-meta.json"
    fi
  fi
done

echo "Done"
