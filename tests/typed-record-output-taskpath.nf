#!/usr/bin/env nextflow
/*
 * Copyright 2013-2026, Seqera Labs
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

nextflow.enable.types = true

workflow {
    sample = record(
        id: 'alpha',
        file: file("${projectDir}/typed-record-output-taskpath/input.txt")
    )

    result = P2(P1(channel.of(sample)))
    result.view { it -> it.text.trim() }
}

process P1 {
    input:
    sample: SampleRecord

    output:
    sample + record(file2: file('test2.txt'))

    script:
    """
    echo 'test2' > test2.txt
    """
}

process P2 {
    input:
    sample: SampleRecord

    output:
    file('combined.txt')

    script:
    """
    test '${sample.file}' = 'input.txt'
    test '${sample.file2}' = 'test2.txt'

    cat '${sample.file}' '${sample.file2}' > combined.txt
    """
}

record SampleRecord {
    id: String
    file: Path
    file2: Path?
}
