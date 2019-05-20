#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.db = "$baseDir/blast-db/tiny"
params.query = "$baseDir/data/sample.fa"
params.chunk = 1 

db = file(params.db)

/* 
 * Extends a BLAST query for each entry in the 'chunks' channel 
 */
process blast {
    input:
    file 'query.fa'

    output:
    file top_hits

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
    file top_hits

    output:
    file sequences

    "blastdbcmd -db ${db} -entry_batch top_hits | head -n 10 > sequences"
}


/*
 * Aligns a T-Coffee MSA and print it 
 */
process align {
    echo true

    input:
    file all_seq

    "t_coffee $all_seq 2>/dev/null | tee align_result"
}

/*
 * main flow
 */
Channel.fromPath(params.query) |
        splitFasta(by: params.chunk) |
        blast |
        extract |
        collectFile(name:'all_seq') |  // Collect all hits to a single file called  'all_seq'
        align

