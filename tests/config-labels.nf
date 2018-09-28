#!/usr/bin/env nextflow
/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

echo true

process alpha {
    /
    echo alpha memry: ${task.memory}
    echo alpha queue: ${task.queue}
    /
}

process beta {
    label 'small'

    /
    echo beta memry: ${task.memory}
    echo beta queue: ${task.queue}
    /
}

process delta {
    label 'big'

    /
    echo delta memry: ${task.memory}
    echo delta queue: ${task.queue}
    /
}

process gamma {
    label 'big'
    memory 4.GB
    queue 'foo'

    /
    echo gamma memry: ${task.memory}
    echo gamma queue: ${task.queue}
    /
}
