#!/usr/bin/env nextflow

params.db = "$HOME/tools/blast-db/pdb/pdb"
params.query = "$HOME/sample.fa"
params.chunkSize = 1

DB=params.db
seq = new Channel()

inputFile = new File(params.query)
inputFile.chunkFasta( params.chunkSize ) { seq << it }

task {
    input '-': seq
    output blastResult

    """
    cat - | blastp -db $DB -query - -outfmt 6 > blastResult
    """
}

merge {
    input blastResult
    output allBlast

    """
    cat ${blastResult} >> allBlast
    """
}

task {
    input allBlast

    """
    sort $allBlast
    """
}


println read(result)
