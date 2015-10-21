#!/bin/bash
TARGET=../../nextflow-website/assets
LATEST=$TARGET/docs/latest/

mkdir -p $TARGET/docs/latest/

rsync -r _build/html/* $LATEST

( cd ../../nextflow-website; jbake)
