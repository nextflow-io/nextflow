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

nextflow.enable.types = true

params {
    input: Path
    prefix: String = 'my'
    names: List<String> = ['alpha', 'delta', 'gamma', 'omega']
}

process foo {
    stageInMode 'copy'

    input:
    tuple(name: String, txt: Path)
    prefix: String

    stage:
    stageAs txt, "${prefix}_${name}.txt"

    output:
    file("${prefix}_${name}.txt")

    script:
    """
    echo World >> ${prefix}_${name}.txt
    """
}

workflow {
    ch_names = channel.fromList(params.names)
    ch_inputs = ch_names.combine( channel.of(params.input) )
    ch_foo = foo( ch_inputs, params.prefix )
    ch_foo.subscribe { it ->
        println "~ Saving ${it.name}"
        it.copyTo('.')
    }
}
