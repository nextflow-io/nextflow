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

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TimeResourceTest extends Specification {

    def 'should create a memory resource' () {

        final _1_h = '1 h' as Duration
        final _2_h = '2 h' as Duration
        final _4_h = '4 h' as Duration

        when:
        def time = new TimeResource(VALUE)
        then:
        time.request == REQ
        time.limit == LIM

        where:
        VALUE                         | REQ      | LIM
        _1_h                          | _1_h     | _1_h
        [request: _2_h]               | _2_h     | null
        [limit: _4_h]                 | _4_h     | _4_h
        [request: _2_h, limit: _4_h]  | _2_h     | _4_h
    }
}
