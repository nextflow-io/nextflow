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

package nextflow.util

import spock.lang.Specification

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class MemoryUnitTest extends Specification {

    def 'should convert mem unit to string'() {

        expect:
        new MemoryUnit(1000).toString() == '1000 B'
        new MemoryUnit(1024).toString() == '1 KB'
        new MemoryUnit(1500).toString() == '1.5 KB'

        new MemoryUnit(2 * 1024 * 1024).toString() == '2 MB'
        new MemoryUnit(3 * 1024L * 1024L * 1024L).toString() == '3 GB'
        new MemoryUnit(4 * 1024L * 1024L * 1024L * 1024L).toString() == '4 TB'
        new MemoryUnit(5 * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '5 PB'
        new MemoryUnit(6 * 1024L * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '6 EB'

    }


    def 'should create mem unit from string' () {

        expect:
        new MemoryUnit('1').toBytes() == 1
        new MemoryUnit('2B').toBytes() == 2
        new MemoryUnit('3 B').toBytes() == 3

        new MemoryUnit('1KB').toBytes() == 1024
        new MemoryUnit('2 KB').toBytes() == 2 * 1024
        new MemoryUnit('3.5 KB').toBytes() == 3.5 * 1024
        new MemoryUnit('5K').toBytes() == 5 * 1024

        new MemoryUnit('1MB').toBytes() == 1024 * 1024
        new MemoryUnit('2 MB').toBytes() == 2 * 1024 * 1024
        new MemoryUnit('3.5 MB').toBytes() == 3.5 * 1024 * 1024
        new MemoryUnit('3.5 M').toBytes() == 3.5 * 1024 * 1024

        new MemoryUnit('1GB').toBytes() == 1024 * 1024 * 1024L
        new MemoryUnit('2 GB').toBytes() == 2 * 1024 * 1024 * 1024L
        new MemoryUnit('3.5 GB').toBytes() == 3.5 * 1024 * 1024 * 1024L
        new MemoryUnit('4G').toBytes() == 4 * 1024 * 1024 * 1024L

        new MemoryUnit('1TB').toBytes() == 1024 * 1024 * 1024L* 1024L
        new MemoryUnit('2 TB').toBytes() == 2 * 1024 * 1024 * 1024L* 1024L
        new MemoryUnit('3.5 TB').toBytes() == 3.5 * 1024 * 1024 * 1024L* 1024L
        new MemoryUnit('25 TB').toBytes() == 25 * 1024 * 1024 * 1024L* 1024L

        new MemoryUnit('1PB').toBytes() == 1024 * 1024 * 1024L* 1024L* 1024L
        new MemoryUnit('2 PB').toBytes() == 2 * 1024 * 1024 * 1024L* 1024L* 1024L
        new MemoryUnit('3.5 PB').toBytes() == 3.5 * 1024 * 1024 * 1024L* 1024L* 1024L
        new MemoryUnit('35 P').toBytes() == 35 * 1024 * 1024 * 1024L* 1024L* 1024L

        new MemoryUnit('1000 KB').toBytes() == 1000 * 1024

        when:
        new MemoryUnit('1,000 KB')
        then:
        thrown(IllegalArgumentException)
    }

    def 'test getters' () {
        
        expect:
        new MemoryUnit('3.5 PB').bytes == 3.5 * 1024 * 1024 * 1024L* 1024L* 1024L
        new MemoryUnit('3.5 PB').kilo == 3.5 * 1024 * 1024 * 1024L* 1024L
        new MemoryUnit('3.5 PB').mega == 3.5 * 1024 * 1024 * 1024L
        new MemoryUnit('3.5 PB').giga == 3.5 * 1024 * 1024

    }

    def 'test equals and compare' () {

        expect:
        new MemoryUnit('1GB') == new MemoryUnit('1GB')
        new MemoryUnit('1M') < new MemoryUnit('1GB')
        new MemoryUnit('1G') > new MemoryUnit('1M')

    }

    def 'test conversion' () {

        def mem

        when:
        mem = new MemoryUnit('100 M')
        then:
        mem.toGiga() == 0
        mem.toMega() == 100
        mem.toKilo() == 100 * 1024L
        mem.toBytes() == 100 * 1024L * 1024L

        when:
        mem = new MemoryUnit('5G')
        then:
        mem.toGiga() == 5
        mem.toMega() == 5 * 1024L
        mem.toKilo() == 5 * 1024L * 1024L
        mem.toBytes() == 5 * 1024L * 1024L * 1024L


        when:
        mem = new MemoryUnit(100_000)
        then:
        mem.toBytes() == 100_000
        mem.toKilo() == 97  // note: this is floor rounded  (97,65625)
        mem.toMega() == 0
        mem.toGiga() == 0
    }

    def 'should multiply memory'()  {

        expect:
        new MemoryUnit('2GB') * 3 == new MemoryUnit('6GB')
        new MemoryUnit('2GB') * 1.5 == new MemoryUnit('3GB')
        // `multiply` a number by a MemoryUnit is implemented by `NumberDelegatingMetaClass`

    }

    def 'should divide memory'()  {

        expect:
        new MemoryUnit('4GB') / 2 == new MemoryUnit('2GB')
        new MemoryUnit('3GB') / 1.5 == new MemoryUnit('2GB')

    }

    def 'should add memory' () {
        expect:
        new MemoryUnit('1GB') + new MemoryUnit('2GB') == new MemoryUnit('3GB')
    }

    def 'should subtract memory' () {
        expect:
        new MemoryUnit('5GB') - new MemoryUnit('2GB') == new MemoryUnit('3GB')
    }

    def 'should validate groovy truth' () {
        expect:
        !new MemoryUnit(0)
        new MemoryUnit(1)
    }

    def 'should validate to unit method' () {
        expect:
        MemoryUnit.of(STR).toUnit(UNIT) == EXPECT
        
        where:
        STR         | UNIT  | EXPECT
        '2 MB'      | 'B'   | 2 * 1024 * 1024
        '2 MB'      | 'KB'  | 2 * 1024
        '2 MB'      | 'MB'  | 2
        '2 MB'      | 'GB'  | 0
        '3.5 GB'    | 'KB'  | 3.5 * 1024 * 1024
    }

}