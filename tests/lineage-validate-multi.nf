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
 * Lineage validation test for multi-process workflows
 * Verifies validation works with complex pipelines having multiple outputs
 */

params.samples = 'sample1,sample2,sample3'

process GENERATE {
    tag "${sample}"
    
    input:
    val sample

    output:
    tuple val(sample), path("${sample}.txt")

    """
    echo "Data for ${sample}" > ${sample}.txt
    echo "Line 2" >> ${sample}.txt
    echo "Line 3" >> ${sample}.txt
    """
}

process ANALYZE {
    tag "${sample}"
    
    input:
    tuple val(sample), path(input_file)

    output:
    tuple val(sample), path("${sample}.analysis.txt")

    """
    wc -l ${input_file} > ${sample}.analysis.txt
    wc -c ${input_file} >> ${sample}.analysis.txt
    """
}

process SUMMARIZE {
    input:
    path analysis_files

    output:
    path 'summary.txt'

    """
    cat ${analysis_files} | sort > summary.txt
    """
}

workflow {
    samples_ch = Channel.of(params.samples.split(','))
    
    GENERATE(samples_ch)
    ANALYZE(GENERATE.out)
    SUMMARIZE(ANALYZE.out.map { it[1] }.collect())
}
