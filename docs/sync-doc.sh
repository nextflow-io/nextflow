#!/bin/bash
TARGET=../../nextflow-web/assets
LATEST=$TARGET/docs/latest/
VERSION=$TARGET/docs/v0.9.1/

mkdir -p $TARGET/docs/latest/
mkdir -p $TARGET/docs/v0.9.1/

rsync -r _build/html/* $LATEST
rsync -r _build/html/* $VERSION

( cd ../../nextflow-web; jbake)
