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
 * Lineage validation test for resume mode
 * Verifies that resumed runs produce equivalent lineage to fresh runs
 */

params.input_value = 'test'

process STEP1 {
    input:
    val x

    output:
    path 'step1.txt'

    """
    echo "Step 1: ${x}" > step1.txt
    sleep 1
    """
}

process STEP2 {
    input:
    path input_file

    output:
    path 'step2.txt'

    """
    cat ${input_file} > step2.txt
    echo "Step 2 complete" >> step2.txt
    sleep 1
    """
}

process STEP3 {
    input:
    path input_file

    output:
    path 'final.txt'

    """
    cat ${input_file} > final.txt
    echo "Step 3 final" >> final.txt
    """
}

workflow {
    STEP1(params.input_value)
    STEP2(STEP1.out)
    STEP3(STEP2.out)
}
