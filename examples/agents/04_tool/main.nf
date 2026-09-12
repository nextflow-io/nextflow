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

// A canonical script task that can run through any Nextflow executor.
process uppercase {
    container 'ubuntu:24.04'
    input:
    text: String

    output:
    result: String = stdout()

    script:
    """
    printf '%s' '${text}' | tr '[:lower:]' '[:upper:]'
    """
}

agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'Call the `uppercase` tool, then reply with only its result.'
    tools 'nf:module_run'

    input:
    request: String

    output:
    answer: String

    prompt:
    """
    ${request}
    """
}

workflow {
    shouty(channel.of('uppercase the word hello'))
        .view { answer -> "ANSWER=${answer}" }
}
