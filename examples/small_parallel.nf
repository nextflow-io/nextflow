params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"

fasta = channel()
queryFile = file(params.query)
queryFile.chunkFasta {
    fasta << it
}

task ('blast') {
    input file:'query.fa', from: fasta
    output top_hits

    """
    blastp -db ${params.db} -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


task ('extract') {
    input top_hits
    output sequences

    "blastdbcmd -db ${params.db} -entry_batch $top_hits > sequences"
}

merge {
    input sequences
    output all_seq

    "cat $sequences > all_seq"
}

task ('align') {
    echo true
    input all_seq

    "t_coffee $all_seq 2>&- | tee align_result"
}
