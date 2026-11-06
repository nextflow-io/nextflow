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

// E1 processes: one native task per item of a channel of n items

params {
    n: Integer = 10000
}

process LIGHT {
    input:
    val x

    output:
    val y

    exec:
    y = x * 2
}

workflow {
    main:
    def t0 = System.nanoTime()
    LIGHT(channel.fromList(1..params.n))
        .reduce(0L) { acc, y -> acc + y }
        .view { v -> "result=${v} dataflow_ms=${(System.nanoTime() - t0).intdiv(1000000)}" }
}
