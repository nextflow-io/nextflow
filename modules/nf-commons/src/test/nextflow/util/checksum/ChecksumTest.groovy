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
package nextflow.util.checksum

import spock.lang.Specification

class ChecksumTest extends Specification {

    def 'should expose algo + value and a canonical hash contribution'() {
        given:
        def ck = new Checksum('crc64nvme', 'abc123')

        expect:
        ck.algo == 'crc64nvme'
        ck.value == 'abc123'
        ck.toHashContribution() == 'crc64nvme:abc123'
    }

    def 'should reject null or blank fields'() {
        when:
        new Checksum(algo, value)

        then:
        thrown(IllegalArgumentException)

        where:
        algo        | value
        null        | 'abc'
        ''          | 'abc'
        'crc64nvme' | null
        'crc64nvme' | ''
    }

    def 'should implement equals and hashCode based on algo + value'() {
        expect:
        new Checksum('crc64nvme', 'abc') == new Checksum('crc64nvme', 'abc')
        new Checksum('crc64nvme', 'abc').hashCode() == new Checksum('crc64nvme', 'abc').hashCode()
        new Checksum('crc64nvme', 'abc') != new Checksum('sha256', 'abc')
        new Checksum('crc64nvme', 'abc') != new Checksum('crc64nvme', 'xyz')
    }
}
