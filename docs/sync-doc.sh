#!/bin/bash
TARGET=../../nextflow-website/assets
if grep -E "^release = '.*edge|.*SNAPSHOT'$" -c conf.py >/dev/null; then
MODE=edge
else
MODE=latest
fi
LATEST=$TARGET/docs/$MODE/

mkdir -p $LATEST

rsync -r _build/html/ $LATEST

( cd ../../nextflow-website; ./jbake)
