#!/usr/bin/env nextflow

params.db = "$HOME/tools/blast-db/pdb/pdb"
params.query = "$HOME/sample.fa"
params.chunkSize = 1

DB = file(params.db)
seq = Channel.create()

inputFile = file(params.query)
inputFile.chunkFasta( params.chunkSize ) { seq << it }

process blast {
    input:
    file 'seq.fa' from seq

    output:
    file 'out' into blast_result

    """
    blastp -db $DB -query seq.fa -outfmt 6 > out
    """
}

process collectResult {
    merge true

    input:
    file blast_result

    output:
    file 'blast_all'

    """
    cat ${blast_result} >> blast_all
    """
}


process sort {
    input:
    file blast_all

    output:
    stdout result

    """
    sort ${blast_all}
    """
}


println result.val
