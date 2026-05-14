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
    result = P1(MAKE_SAMPLE())
    result.view { r -> "${r.file.text.trim()} ${r.file2.text.trim()}" }
}

process MAKE_SAMPLE {
    output:
    record(
        id: 'alpha',
        file: file('input.txt')
    )

    script:
    """
    echo 'test1' > input.txt
    """
}

process P1 {
    input:
    sample: Sample

    output:
    sample + record(file2: file('test2.txt'))

    script:
    """
    echo 'test2' > test2.txt
    """
}

record Sample {
    id: String
    file: Path
    file2: Path?
}
