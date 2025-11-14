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

process matmul1 {
    accelerator 1
    debug true

    input:
    val id

    output:
    val id

    script:
    '''
    echo "matmul1: gpu $CUDA_VISIBLE_DEVICES"
    sleep 5
    '''
}

process matmul2 {
    accelerator 1
    debug true

    input:
    val id

    script:
    '''
    echo "matmul2: gpu $CUDA_VISIBLE_DEVICES"
    sleep 5
    '''
}

workflow {
    matmul1( channel.of(1..8) ) | matmul2
}
