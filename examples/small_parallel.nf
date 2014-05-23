#!/usr/bin/env nextflow

params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"

db = file(params.db)
fasta = Channel.fromPath(params.query).splitFasta()

process blast {
    input:
    file 'query.fa' from fasta

    output:
    file 'top_hits'

    """
    blastp -db ${db} -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


process extract {
    input:
    file top_hits

    output:
    file 'sequences'

    "blastdbcmd -db ${db} -entry_batch top_hits | head -n 10 > sequences"
}


all_seq = sequences.collectFile(name:'all_seq')

process align {
    echo true

    input:
    file all_seq

    "t_coffee $all_seq 2>&- | tee align_result"
}

