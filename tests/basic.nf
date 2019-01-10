#!/usr/bin/env nextflow
/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
 
 
/* 
 * Command line input parameter 
 */
params.in = "$baseDir/data/sample.fa"

/* 
 * Define the input file 
 */
sequences = file(params.in)


/* 
 * split a fasta file in multiple files 
 */
process splitSequences {

    input:
    file 'input.fa' from sequences

    output:
    file 'seq_*' into records

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """

}

/* 
 * Simple reverse the sequences 
 */
process reverse {

    input:
    file x from records
    
    output:
    stdout result

    """
    cat $x | rev
    """
}

/* 
 * print the channel content 
 */
result.subscribe { println it }
