#!/usr/bin/env nextflow

params.db = "$baseDir/blast-db/tiny"
params.query = "$baseDir/data/sample.fa"
params.chunk = 1 

/* 
 * Extends a BLAST query for each entry in the 'chunks' channel 
 */
process blast {
    input:
    path 'query.fa'
    val db

    output:
    path 'top_hits'

    script:
    """
    blastp -db ${db} -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}

/*
 * Find out the top 10 matches returned by the BLAST query
 */ 
process extract {
    input:
    path 'top_hits'
    path db

    output:
    path 'sequences'

    script:
    "blastdbcmd -db ${db} -entry_batch top_hits | head -n 10 > sequences"
}


/*
 * Aligns a T-Coffee MSA and print it 
 */
process align {
    debug true

    input:
    path all_seq

    script:
    "t_coffee $all_seq 2>/dev/null | tee align_result"
}

/*
 * main flow
 */
workflow {
    ch_fasta = Channel.fromPath(params.query)
        | splitFasta(by: params.chunk, file:true)

    ch_sequences = blast(ch_fasta, params.db)

    extract(ch_sequences, params.db)
        | collectFile(name:'all_seq') // Collect all hits to a single file called  'all_seq'
        | align
}
