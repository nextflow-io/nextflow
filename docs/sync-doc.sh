#!/bin/bash
TARGET=../../nextflow-website/assets
if egrep "^release = '.*edge'$" -c conf.py >/dev/null; then
MODE=edge
else 
MODE=latest
fi
LATEST=$TARGET/docs/$MODE/

mkdir -p $TLATEST

rsync -r _build/html/* $LATEST

( cd ../../nextflow-website; jbake)
