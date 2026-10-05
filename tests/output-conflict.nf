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

params.conflict = false

process MAKE {
    input:
    val id

    output:
    tuple val(id), path('versions.yml'), path('a.txt'), path('b.txt')

    script:
    """
    echo "${id}: 1.0" > versions.yml
    echo "a ${id}" > a.txt
    echo "b ${id}" > b.txt
    """
}

workflow {
    main:
    ch = MAKE(channel.of('x', 'y', 'z'))
        .map { id, versions, a, b -> [id: id, versions: versions, a: a, b: b] }

    publish:
    versions = ch.map { r -> r.versions }
    reports = ch
}

output {
    // every value publishes a file to the same target, which is allowed
    versions {
        path 'info'
        index {
            path 'versions.csv'
        }
    }

    // the same value publishes two files to the same target when `--conflict` is set
    reports {
        path { r ->
            r.a >> "reports/${r.id}.txt"
            r.b >> (params.conflict ? "reports/${r.id}.txt" : "reports/${r.id}.b.txt")
        }
    }
}
