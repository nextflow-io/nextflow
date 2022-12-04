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
workflow {
    alpha()
    beta()
    delta()
    gamma()
}

process alpha {
    debug true
    /
    echo alpha memry: ${task.memory}
    echo alpha queue: ${task.queue}
    /
}

process beta {
    debug true
    label 'small'

    /
    echo beta memry: ${task.memory}
    echo beta queue: ${task.queue}
    /
}

process delta {
    debug true
    label 'big'

    /
    echo delta memry: ${task.memory}
    echo delta queue: ${task.queue}
    /
}

process gamma {
    debug true
    label 'big'
    memory 40.MB
    queue 'foo'

    /
    echo gamma memry: ${task.memory}
    echo gamma queue: ${task.queue}
    /
}
