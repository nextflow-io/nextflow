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

package nextflow.executor.res

import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MemoryResourceTest extends Specification {

    def 'should create a memory resource' () {

        final _1_GB = '1 GB' as MemoryUnit
        final _2_GB = '2 GB' as MemoryUnit
        final _4_GB = '4 GB' as MemoryUnit

        when:
        def mem = new MemoryResource(VALUE)
        then:
        mem.request == REQ
        mem.limit == LIM

        where:
        VALUE                           | REQ       | LIM
        _1_GB                           | _1_GB     | _1_GB
        [request: _2_GB]                | _2_GB     | null
        [limit: _4_GB]                  | _4_GB     | _4_GB
        [request: _2_GB, limit: _4_GB]  | _2_GB     | _4_GB
    }
}
