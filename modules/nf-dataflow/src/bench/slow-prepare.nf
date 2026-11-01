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

// E3 slow preparation: one native task per item whose `when` guard blocks for `ms` milliseconds,
// standing in for slow staging, hashing or cache lookups, next to a pure dataflow pipeline

params {
    n: Integer = 1000
    ms: Integer = 20
    side: Integer = 100000
}

process SLOW_PREPARE {
    input:
    val x

    output:
    val y

    when:
    { -> Thread.sleep(params.ms); true }()

    exec:
    y = x * 2
}

workflow {
    main:
    def t0 = System.nanoTime()
    SLOW_PREPARE(channel.fromList(1..params.n))
        .reduce(0L) { acc, y -> acc + y }
        .view { v -> "result=${v} dataflow_ms=${(System.nanoTime() - t0).intdiv(1000000)}" }

    channel.fromList(1..params.side)
        .map { x -> x * 2 }
        .reduce(0L) { acc, y -> acc + y }
        .view { v -> "side=${v} side_ms=${(System.nanoTime() - t0).intdiv(1000000)}" }
}
