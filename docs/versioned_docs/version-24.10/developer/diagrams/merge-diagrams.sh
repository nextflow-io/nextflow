#!/bin/bash

packages=()
packages+=("nextflow")
# packages+=("nextflow.ast")
packages+=("nextflow.cache")
packages+=("nextflow.cli")
# packages+=("nextflow.cloud.aws")
# packages+=("nextflow.cloud.aws.nio")
# packages+=("nextflow.cloud.azure")
# packages+=("nextflow.cloud.google")
packages+=("nextflow.config")
# packages+=("nextflow.container")
packages+=("nextflow.dag")
# packages+=("nextflow.executor")
# packages+=("nextflow.extension")
# packages+=("nextflow.ga4gh")
# packages+=("nextflow.k8s")
# packages+=("nextflow.plugin")
packages+=("nextflow.processor")
# packages+=("nextflow.scm")
packages+=("nextflow.script")
# packages+=("nextflow.secret")
# packages+=("nextflow.trace")

outfile="nextflow-merged.mmd"

echo "classDiagram" > ${outfile}

for package in "${packages[@]}"; do
    echo "${package}"

    tail -n +2 "${package}.mmd" >> ${outfile}
    echo >> ${outfile}
done
