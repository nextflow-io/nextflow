/*
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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AcceleratorResourceTest extends Specification {

    def 'should create a gpu resource' () {

        when:
        def acc = new AcceleratorResource(VALUE)
        then:
        acc.type == TYPE
        acc.request == REQ
        acc.limit == LIM
        acc.type == TYPE
        acc.runtime == RUNTIME

        where:
        VALUE                       | REQ   | LIM   | TYPE  | RUNTIME
        1                           | 1     | 1     | null  | null
        5                           | 5     | 5     | null  | null
        [request: 2]                | 2     | null  | null  | null
        [limit: 4]                  | 4     | 4     | null  | null
        [request: 2, limit: 4]      | 2     | 4     | null  | null
        [request: 2, limit: 4]      | 2     | 4     | null  | null
        [limit: 3, type: 'nvidia']  | 3     | 3     | 'nvidia' | null
        [limit: 3, runtime: 'foo']  | 3     | 3     | null  | 'foo'
    }
}
