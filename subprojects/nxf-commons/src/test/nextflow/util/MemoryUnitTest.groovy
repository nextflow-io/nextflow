/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util

import spock.lang.Specification

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class MemoryUnitTest extends Specification {

    def 'test toString'() {

        expect:
        new MemoryUnit(1024).toString() == '1 KB'
        new MemoryUnit(2 * 1024 * 1024).toString() == '2 MB'
        new MemoryUnit(3 * 1024L * 1024L * 1024L).toString() == '3 GB'
        new MemoryUnit(4 * 1024L * 1024L * 1024L * 1024L).toString() == '4 TB'
        new MemoryUnit(5 * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '5 PB'
        new MemoryUnit(6 * 1024L * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '6 EB'

    }


    def 'test create' () {

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

}