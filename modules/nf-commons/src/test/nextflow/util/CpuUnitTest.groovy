/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.util


import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CpuUnitTest extends Specification {

    @Unroll
    def 'should parse cpu string'() {
        expect:
        CpuUnit.parseCpuMillis(VALUE) == MILLIS

        where:
        VALUE   | MILLIS
        '1'     | 1_000
        '8'     | 8_000
        '256'   | 256_000
        and:
        '10m'   | 10
        '999m'  | 999
        '1100m' | 1100
        and:
        '1.0'    | 1_000
        '1.1'    | 1_100
        '8.0'    | 8_000
    }

    @Unroll
    def 'should parse cpu number'() {
        expect:
        CpuUnit.parseCpuMillis(VALUE) == MILLIS

        where:
        VALUE   | MILLIS
        1       | 1_000
        8       | 8_000
        and:
        1.0     | 1_000
        1.1     | 1_100
        8.0     | 8_000
    }

    def 'should validate equals and hash code' () {
        given:
        def c1 = CpuUnit.of(100)
        def c2 = CpuUnit.of(100)
        def c3 = CpuUnit.of(300)

        expect:
        c1 == c2 
        c1 != c3
        and:
        c1.hashCode() == c2.hashCode()
        c1.hashCode() != c3.hashCode()

    }

    def 'should get cores' () {
        expect:
        CpuUnit.of(VALUE).toCores() == CORES

        where:
        VALUE       | CORES
        1           | 1
        2           | 2
        16          | 16
        256         | 256
        and:
        '100m'      | 1
        '200m'      | 1
        '999m'      | 1
        '1000m'     | 1
        and:
        0.1         | 1
        0.99        | 1
        1.1         | 2
    }

    def 'should get milils' () {
        expect:
        CpuUnit.of(VALUE).toMillis() == MILLIS

        where:
        VALUE       | MILLIS
        1           | 1_000
        2           | 2_000
        16          | 16_000
        256         | 256_000
        and:
        '100m'      | 100
        '200m'      | 200
        '999m'      | 999
        '1000m'     | 1000
        and:
        0.1         | 100
        0.99        | 990
        1.1         | 1_100
    }

    def 'should get decimal string' () {
        expect:
        CpuUnit.of(VALUE).toDecimalString() == STRING

        where:
        VALUE       | STRING
        1           | '1.0'
        2           | '2.0'
        16          | '16.0'
        256         | '256.0'
        2100        | '2100.0'
        and:
        '100m'      | '0.1'
        '200m'      | '0.2'
        '999m'      | '1.0'
        '1000m'     | '1.0'
        and:
        0.01        | '0.1'
        0.1         | '0.1'
        0.9         | '0.9'
        0.91        | '1.0'
        0.99        | '1.0'
        1.1         | '1.1'
    }

    def 'test equals and compare' () {

        expect:
        CpuUnit.of(4) == CpuUnit.of('4000m')
        CpuUnit.of(1.2) < CpuUnit.of(10)
        CpuUnit.of('5100m') > CpuUnit.of(3.6)

    }
}
