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

params.in = "$baseDir/data/sample.fa"
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process split {
    input:
    file 'query.fa' from file(params.in)

    output:
    file 'seq_*' into splits

    """
    $SPLIT query.fa '%^>%' '/^>/' '{*}' -f seq_
    """
}


process printTwo {
    echo true

    input:
    file 'chunk' from splits

    output:
    file 'chunk1:chunk3' into two_chunks mode flatten

    """
    cat chunk* | rev
    """

}

process printLast {
    echo true

    input:
    file 'chunk' from two_chunks

    output:
    file 'chunk' into result

    """
    cat chunk
    """
}
