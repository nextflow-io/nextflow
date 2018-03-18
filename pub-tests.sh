#!/bin/bash
set -x
#
# Make test results available through
#    https://www.nextflow.io/tests/travis/index.html
#
export AWS_ACCESS_KEY_ID=$NXF_AWS_ACCESS
export AWS_SECRET_ACCESS_KEY=$NXF_AWS_SECRET
export AWS_DEFAULT_REGION=eu-west-1

aws --region $AWS_DEFAULT_REGION s3 rm --quiet --recursive s3://www.nextflow.io/tests/$1

for x in $(find . -path \*build/reports/tests/test); do
aws --region $AWS_DEFAULT_REGION s3 sync $x s3://www.nextflow.io/tests/$1/${x#./} \
 --cache-control max-age=0 \
 --metadata-directive REPLACE \
 --storage-class REDUCED_REDUNDANCY \
 --acl public-read \
 --delete
done