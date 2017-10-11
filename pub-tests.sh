#!/usr/bin/env bash

#
# Make test results available through
#    https://www.nextflow.io/tests/travis/index.html
#
export AWS_ACCESS_KEY_ID=$NXF_AWS_ACCESS
export AWS_SECRET_ACCESS_KEY=$NXF_AWS_SECRET
export AWS_DEFAULT_REGION=eu-west-1

aws --region eu-west-1 s3 sync $1 s3://www.nextflow.io/tests/$2/ \
 --cache-control max-age=0 \
 --metadata-directive REPLACE \
 --storage-class REDUCED_REDUNDANCY \
 --acl public-read \
 --delete
