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

seqs = Channel.fromPath("$baseDir/data/p{1,2,3}.fa")

process foo {
    echo true

    input:
    file 'dir1/link_*.fasta' from seqs.toList()

    output:
    file 'dir2/*' into result mode flatten

    '''
    ls dir1 | sort
    mkdir dir2
    echo hello > dir2/alpha.txt
    echo world > dir2/beta.txt
    '''
}

result.println { it.name }
