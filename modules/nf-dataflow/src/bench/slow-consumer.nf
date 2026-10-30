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

// E2 slow consumer: n records into a `maxForks 1` process with a tiny native task
// (port of the `mem` workload of nf prototypes/0014-dataflow). Without backpressure
// the source binds every record up front and the process turns each one into a
// pending task, so memory grows with n.

params {
    n: Integer = 1000000
}

process SLOW {
    maxForks 1

    input:
    val record

    output:
    val y

    exec:
    y = record[1]
}

workflow {
    main:
    def t0 = System.nanoTime()
    def records = channel.fromList(1..params.n)
        .map { i -> ["sample_${i}", i as long, "/data/reads/sample_${i}.fq.gz"] }
    SLOW(records)
        .count()
        .view { v -> "result=${v} dataflow_ms=${(System.nanoTime() - t0).intdiv(1000000)}" }
}
