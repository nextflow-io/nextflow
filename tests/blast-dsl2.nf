#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.db = "$baseDir/blast-db/tiny"
params.query = "$baseDir/data/sample.fa"
params.chunkSize = 1

DB = file(params.db)

process blast {
    input:
    file 'seq.fa'

    output:
    file 'out'

    """
    blastp -db $DB -query seq.fa -outfmt 6 > out
    """
}

process sort {
    input:
    file 'hits_*'

    output:
    stdout result

    """
    sort hits_*
    """
}

Channel.fromPath(params.query) |
        splitFasta( by: params.chunkSize ) |
        blast |
        collect |
        sort |
        subscribe { println it }
