#!/bin/bash

docker buildx build \
    --push \
    --platform linux/amd64,linux/arm64 \
    --tag nextflow/sphinx:5.3.0 \
    .
