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

nextflow.preview.types = true

include { SampleRecord } from './typed-record-staging/records/main.nf'
include { CHECK_SAMPLE } from './typed-record-staging/modules/check/main.nf'

workflow {
    resource = new SampleRecord(
        id: 'alpha',
        path: file("${projectDir}/typed-record-staging/input.txt")
    )

    result = CHECK_SAMPLE(channel.of(resource))
    result.view { it -> it.text.trim() }
}