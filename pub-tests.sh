#!/bin/bash
set -x
#
# Make test results available through
#    https://www.nextflow.io/tests/travis/index.html
#
export AWS_ACCESS_KEY_ID=${NXF_AWS_ACCESS}
export AWS_SECRET_ACCESS_KEY=${NXF_AWS_SECRET}
export AWS_DEFAULT_REGION=eu-west-1

# use the java version as target context
ver=${TEST_JDK:-8}
base=s3://www2.nextflow.io/tests/$1/jdk$ver

# delete previous run
aws --region $AWS_DEFAULT_REGION s3 rm --only-show-errors --recursive $base

# upload unit test results
for x in $(find . -path \*build/reports/tests/test); do
aws --region $AWS_DEFAULT_REGION s3 sync $x $base/${x#./} \
 --only-show-errors \
 --cache-control max-age=0 \
 --metadata-directive REPLACE \
 --storage-class STANDARD \
 --acl public-read \
 --delete
done

# upload integration test results
aws --region $AWS_DEFAULT_REGION s3 sync --only-show-errors tests/checks $base/integration