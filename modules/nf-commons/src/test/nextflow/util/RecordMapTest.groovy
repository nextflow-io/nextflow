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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class RecordMapTest extends Specification {

    def testRecord() {
        when:
        def record = new RecordMap(id: '1', fastq_1: '1_1.fastq', fastq_2: '1_2.fastq')
        then:
        record.id == '1'
        record.fastq_1 == '1_1.fastq'
        record.fastq_2 == '1_2.fastq'

        // when:
        // record.foo = 'bar'
        // then:
        // thrown(UnsupportedOperationException)

        // when:
        // record.put(foo: 'bar')
        // then:
        // thrown(UnsupportedOperationException)

        when:
        record.remove('id')
        then:
        thrown(UnsupportedOperationException)

        expect:
        record.subMap(['id']) == new RecordMap(id: '1')

    }
}
