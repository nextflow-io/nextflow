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

package nextflow.executor.res

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AcceleratorResourceTest extends Specification {

    def 'should create an accelerator resource' () {

        when:
        def acc = new AcceleratorResource(VALUE)
        then:
        acc.type == TYPE
        acc.request == REQ
        acc.limit == LIM
        acc.type == TYPE
        acc.runtime == RUNTIME
        acc.cuda == CUDA

        where:
        VALUE                       | REQ   | LIM   | TYPE  | RUNTIME | CUDA
        1                           | 1     | 1     | null  | null  | null
        5                           | 5     | 5     | null  | null  | null
        [request: 2]                | 2     | null  | null  | null  | null
        [limit: 4]                  | 4     | 4     | null  | null  | null
        [request: 2, limit: 4]      | 2     | 4     | null  | null  | null
        [limit: 3, type: 'nvidia']  | 3     | 3     | 'nvidia' | null | null
        [limit: 3, runtime: 'foo']  | 3     | 3     | null  | 'foo' | null
        [request: 1, cuda: '12.2']  | 1     | null  | null  | null  | '12.2'
        [request: 1, type: 'nvidia-l4', cuda: '12.10'] | 1 | null | 'nvidia-l4' | null | '12.10'
    }

    def 'should reject an invalid cuda version' () {
        when:
        new AcceleratorResource([request: 1, cuda: VALUE])
        then:
        def err = thrown(IllegalArgumentException)
        err.message.contains('Invalid `accelerator` directive `cuda` option')

        where:
        VALUE << [12.2, 12, '12', '12.x', 'v12.2', '12.2.1', '']
    }

    def 'should throw error for invalid resource' () {
        when:
        def acc = new AcceleratorResource([count: 1, type: 'nvidia-tesla-t4'])
        then:
        def err = thrown(IllegalArgumentException)
        err.message.contains('Invalid `accelerator` directive value -- `request` or `limit` is required')
    }
}
