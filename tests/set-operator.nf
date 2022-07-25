#!/usr/bin/env nextflow
/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
nextflow.enable.dsl=2

Channel
    .of ( 'a', 'b', 'c' )
    .set { log1 }
    .map { it * 2 }
    .set { log2 }
    .map { it.toUpperCase() }
    .view { "Result: $it" }

log1.view { "Log 1: $it" }
log2.view { "Log 2: $it" }
