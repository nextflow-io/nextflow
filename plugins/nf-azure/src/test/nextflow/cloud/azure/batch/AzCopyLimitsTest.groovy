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

package nextflow.cloud.azure.batch

import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Jonathan Manning <jonathan.manning@seqera.io>
 */
class AzCopyLimitsTest extends Specification {

    def 'should not set limits when the task memory is not specified' () {
        expect:
        AzCopyLimits.of(null, 4, 4) == null
    }

    @Unroll
    def 'should size azcopy for #MEM' () {
        when:
        def limits = AzCopyLimits.of(MemoryUnit.of(MEM), MAX_PARALLEL, BLOCK_MB)

        then:
        limits.parallelTransfers == PARALLEL
        limits.bufferMb == BUFFER_MB
        limits.bufferGb == BUFFER_GB

        where:
        MEM      | MAX_PARALLEL | BLOCK_MB | PARALLEL | BUFFER_MB | BUFFER_GB
        '128 MB' | 4            | 4        | 1        | 32        | '0.031'
        '512 MB' | 4            | 4        | 1        | 64        | '0.062'
        '1 GB'   | 4            | 4        | 2        | 128       | '0.125'
        '2 GB'   | 4            | 4        | 4        | 160       | '0.156'
        '4 GB'   | 4            | 4        | 4        | 416       | '0.406'
        '16 GB'  | 4            | 4        | 4        | 1952      | '1.906'
        '16 GB'  | 2            | 4        | 2        | 3968      | '3.875'
        '512 MB' | 4            | 100      | 1        | 400       | '0.39'
    }
}
