params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"

fasta = channel()
queryFile = file(params.query)
queryFile.chunkFasta {
    fasta << it
}

process blast {
    input:
    file 'query.fa' using fasta

    output:
    file top_hits

    """
    blastp -db ${params.db} -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


process extract {
    input:
    val top_hits

    output:
    file sequences

    "blastdbcmd -db ${params.db} -entry_batch $top_hits > sequences"
}

process all(merge:true) {
    input:
    val sequences

    output:
    file all_seq

    "cat $sequences > all_seq"
}

process align {
    echo true

    input:
    val all_seq

    "t_coffee $all_seq 2>&- | tee align_result"
}
