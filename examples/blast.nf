#!/usr/bin/env nextflow

params.db = "$HOME/tools/blast-db/pdb/pdb"
params.query = "$HOME/sample.fa"
params.chunkSize = 1

DB= file(params.db)
seq = channel()

inputFile = file(params.query)
inputFile.chunkFasta( params.chunkSize ) { seq << it }

process blast {
    input:
    file seq as 'seq.fa'

    output:
    file 'out' to blast_result

    """
    blastp -db $DB -query seq.fa -outfmt 6 > out
    """
}

blast_all = blast_result.toList()

process sort {
    input:
    file blast_all as 'hits_*'

    output:
    stdout result

    """
    sort hits_*
    """
}


result.subscribe { println it }
