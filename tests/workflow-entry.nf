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

// Test direct execution of a named workflow with parameter mapping
nextflow.enable.types = true

record Sample {
    id: String
    data: Path
}

workflow SAY_HELLO {
    take:
    samples: Channel<Sample>
    prefix: String
    limit: Integer?

    main:
    greetings = samples.map { s -> "${prefix}, ${s.id} (${s.data.name})" }

    emit:
    greetings: Channel<String> = greetings
    info: Value<String> = "prefix=${prefix} limit=${limit}"
}
