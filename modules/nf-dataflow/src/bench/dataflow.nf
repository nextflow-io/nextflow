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

// E1 dataflow: n items through map, filter, flatMap, map and reduce
// (port of nf bench/run/dataflow.nf)

params {
    n: Integer = 1000000
}

workflow {
    main:
    def t0 = System.nanoTime()
    channel.fromList(0..<params.n)
        .map { x -> x + 1 }
        .filter { x -> x % 3 != 0 }
        .flatMap { x -> [x, x] }
        .map { x -> x * 2 }
        .reduce(0L) { acc, x -> acc + x }
        .view { v -> "result=${v} dataflow_ms=${(System.nanoTime() - t0).intdiv(1000000)}" }
}
