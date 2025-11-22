$NXF_CMD run \
    rnaseq-nf \
    -profile batch,s3-data \
    -with-wave \
    -with-fusion \
    -process.scratch false

