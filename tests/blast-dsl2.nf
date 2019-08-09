#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.db = "$baseDir/blast-db/tiny"
params.query = "$baseDir/data/sample.fa"
params.chunkSize = 1

DB = file(params.db)

process blast {
    input:
    path 'seq.fa'

    output:
    path 'out'

    """
    blastp -db $DB -query seq.fa -outfmt 6 > out
    """
}

process sort {
    input:
    path 'hits_*'

    output:
    stdout()

    """
    sort hits_*
    """
}


workflow {
    Channel.fromPath(params.query) |
            splitFasta( by: params.chunkSize, file:true ) |
            blast |
            collect |
            sort |
            subscribe { println it }
}
