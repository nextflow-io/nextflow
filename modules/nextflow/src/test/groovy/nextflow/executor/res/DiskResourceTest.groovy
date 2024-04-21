/*
 * Copyright 2013-2024, Seqera Labs
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

    static final _100_GB = MemoryUnit.of('100GB')
    static final _375_GB = MemoryUnit.of('375GB')

    def 'should create a disk resource' () {

        when:
        def disk = new DiskResource(VALUE)
        then:
        disk.request == REQ
        disk.type == TYPE

        where:
        VALUE                                  | REQ      | TYPE
        _100_GB                                | _100_GB  | null
        [request: _100_GB]                     | _100_GB  | null
        [request: _375_GB, type: 'local-ssd']  | _375_GB  | 'local-ssd'
    }

    def 'should return a disk resource with the specified request' () {
        expect:
        new DiskResource(request: _100_GB).withRequest(_375_GB) == new DiskResource(request: _375_GB)
        new DiskResource(request: _100_GB, type: 'ssd').withRequest(_375_GB) == new DiskResource(request: _375_GB, type: 'ssd')
    }
}
