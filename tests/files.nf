#!/usr/bin/env nextflow
/*
 * Copyright 2013-2024, Seqera Labs
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
    path 'query.fa'

    output:
    path 'seq_*'

    """
    $SPLIT query.fa '%^>%' '/^>/' '{*}' -f seq_
    """
}


process printTwo {
    debug true

    input:
    path 'chunk'

    output:
    file 'chunk1:chunk3'

    """
    cat chunk* | rev
    """

}

process printLast {
    debug true

    input:
    file 'chunk'

    output:
    file 'chunk'

    """
    cat chunk
    """
}

workflow {
  split(params.in) | printTwo | flatten | printLast
}
