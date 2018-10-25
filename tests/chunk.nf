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

/*
 * Show how get data from the 'stdin' virtual file
 *
 * The data is supposed to be a FASTA file which is splitted
 * is chunks containing a single sequence
 *
 * The number of sequences in each each is controlled by the command
 * line parameter '--chunkSize' (--chunk-size is a synonym for the same)
 */


params.chunkSize = 1
sequences = Channel.create()

stdin.splitFasta( by: params.chunkSize, into: sequences )

process foo {
    echo true

    input:
    stdin sequences

    "cat -"
}
