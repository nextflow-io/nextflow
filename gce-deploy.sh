#!/bin/bash
# TODO: This is a temporary workaround to deploy local nextflow binary go GCE bucket

for build in build/releases/nextflow-*-all; do
  if [ "$build" = 'build/releases/nextflow-*-all' ]; then
    echo "No build found.  run 'make compile pack install'"
    exit 1
  fi
  version=$(echo $build | sed 's/^.*nextflow-//' | sed 's/-all$//')
  gsutil cp $build gs://nextflow-nextcode-dev/releases/v${version}/nextflow
  #gsutil cp $build gs://sekretname/releases/v${version}/nextflow
done

echo "To have nextflow cluster download binaries from google bucket:"
#echo "export NEXTFLOW_DOWNLOAD_URL=http://sekretname.storage.googleapis.com/releases"
echo "export NEXTFLOW_DOWNLOAD_URL=http://nextflow-nextcode-dev.storage.googleapis.com/releases"
echo "before running nextflow cloud create"