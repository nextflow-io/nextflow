#!/usr/bin/env nextflow

params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"

db = file(params.db)
query = file(params.query)

process blast {
    input:
    file query
    
    output:
    file 'top_hits'

    """
    blastp -db $db -query $query -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


process extractTopHits {
    input:
    file top_hits

    output:
    file 'sequences'

    "blastdbcmd -db ${db} -entry_batch $top_hits > sequences"
}

process align {
    echo true

    input:
    file sequences

    "t_coffee $sequences 2>&- | tee align_result"
}
