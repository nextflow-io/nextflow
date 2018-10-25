#!/usr/bin/env nextflow
/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

params.db = "$baseDir/blast-db/tiny"
params.query = "$baseDir/data/sample.fa"
params.chunk = 1 

db = file(params.db)
chunks = Channel
            .fromPath(params.query)
            .splitFasta(by: params.chunk)

/* 
 * Extends a BLAST query for each entry in the 'chunks' channel 
 */
process blast {
    input:
    file 'query.fa' from chunks

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
 * Collect all hits to a single file called  'all_seq'
 */ 
all_seq = sequences.collectFile(name:'all_seq')

/*
 * Aligns a T-Coffee MSA and print it 
 */
process align {
    echo true

    input:
    file all_seq

    "t_coffee $all_seq 2>/dev/null | tee align_result"
}

