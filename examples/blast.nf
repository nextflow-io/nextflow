#!/usr/bin/env nextflow

params.db = "$HOME/tools/blast-db/pdb/pdb"
params.query = "$HOME/sample.fa"
params.chunkSize = 1

DB=params.db
seq = channel()

inputFile = file(params.query)
inputFile.chunkFasta( params.chunkSize ) { seq << it }

task {
    stdin seq
    output file: 'out', into: blast_result

    """
    cat - | blastp -db $DB -query - -outfmt 6 > out
    """
}

merge {
    input blast_result
    output blast_all

    """
    cat ${blast_result} >> blast_all
    """
}

task {
    input blast_all

    """
    sort $blast_all
    """
}


println read(result)
