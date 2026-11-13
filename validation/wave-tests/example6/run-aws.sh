$NXF_CMD run \
    rnaseq-nf \
    -profile batch,s3-data \
    -w s3://nextflow-ci-oss/work \
    -with-wave \
    -with-fusion \
    -process.scratch false
