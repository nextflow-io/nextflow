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
class DiskResourceTest extends Specification {

    static final _1_GB = MemoryUnit.of('1GB')
    static final _2_GB = MemoryUnit.of('2GB')
    static final _4_GB = MemoryUnit.of('4GB')

    def 'should create a disk resource' () {

        when:
        def mem = new DiskResource(VALUE)
        then:
        mem.request == REQ
        mem.limit == LIM

        where:
        VALUE                           | REQ    | LIM
        _1_GB                           | _1_GB  | _1_GB
        [request: _2_GB]                | _2_GB  | null
        [limit: _4_GB]                  | _4_GB  | _4_GB
        [request: _2_GB, limit: _4_GB]  | _2_GB  | _4_GB
    }

}
