$NXF_CMD run \
    rnaseq-nf \
    -profile google-batch,gs-data \
    -with-wave \
    -with-fusion \
    -process.scratch false

