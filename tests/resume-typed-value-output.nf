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

/*
 * A typed (v2) *script* process mixing genuine in-memory VALUE outputs with
 * work-dir derived ones (env/file/stdout). A script task always re-evaluates
 * its body before the cache is consulted, so the task context is rebuilt on
 * resume and the value outputs must be emitted identically without the
 * context being persisted in the cache DB.
 */
process foo {
    input:
    message: String

    output:
    record(
        upper: message.toUpperCase(),
        size: message.length(),
        out_env: env('MESSAGE'),
        out_file: file('message.txt'),
        out_std: stdout()
    )

    script:
    """
    export MESSAGE='${message}-env'
    echo '${message}' > message.txt
    printf 'out:${message}'
    """
}

/*
 * A typed (v2) *exec* process: the body IS the task execution, so on a cache hit it is
 * never re-run and its value output can only come from the persisted task context.
 */
process bar {
    input:
    message: String

    output:
    rev

    exec:
    rev = message.reverse()
}

workflow {
    reversed = bar(channel.of('world'))
    reversed.view { v -> "rev=${v}" }

    res = foo(channel.of('hello'))
    res.view { r -> "upper=${r.upper}" }
    res.view { r -> "size=${r.size}" }
    res.view { r -> "env=${r.out_env}" }
    res.view { r -> "file=${r.out_file.text.trim()}" }
    res.view { r -> "std=${r.out_std}" }
}
