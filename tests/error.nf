#!/usr/bin/env nextflow
/*
 * Copyright 2013-2023, Seqera Labs
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

process task1 {
    maxForks 4
    errorStrategy 'ignore'

    input:
    val x

    script:
    "echo $x; exit 1"
}

process task2 {
    maxForks 4

    input:
    val x

    script:
    "echo $x"
 }

workflow {
  channel.of(1,2,3) | task1
  channel.of(4,5,6) | task2
}
