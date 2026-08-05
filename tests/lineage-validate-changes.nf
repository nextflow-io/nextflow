#!/usr/bin/env nextflow
/*
 * Copyright 2013-2025, Seqera Labs
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
 * Lineage validation test for detecting changes
 * Verifies that validation detects parameter and output differences
 */

params.input_text = 'default text'
params.multiplier = 1

process TRANSFORM {
    input:
    val text
    val multiplier

    output:
    path 'output.txt'

    """
    # Repeat text based on multiplier
    for i in \$(seq 1 ${multiplier}); do
        echo "${text}"
    done > output.txt
    """
}

process CHECKSUM {
    input:
    path input_file

    output:
    path 'checksum.txt'

    """
    md5sum ${input_file} | cut -d' ' -f1 > checksum.txt
    """
}

workflow {
    TRANSFORM(params.input_text, params.multiplier)
    CHECKSUM(TRANSFORM.out)
}
