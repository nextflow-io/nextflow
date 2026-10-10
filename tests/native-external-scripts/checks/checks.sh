#!/usr/bin/env bash
set -euo pipefail

# Native external scripts are an executable contract for
# adr/20260522-native-external-scripts.md. This check is kept with the fixture
# until the proposed script-file syntax and runtime context injection are
# implemented.
: "${NXF_RUN:=nextflow -q run ../main.nf}"

$NXF_RUN | tee .stdout

diff checks/expected.txt .stdout
