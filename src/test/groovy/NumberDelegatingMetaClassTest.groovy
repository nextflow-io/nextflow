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

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NumberDelegatingMetaClassTest extends Specification {


    def 'should parse duration units' () {
        expect:
        1.ms  == Duration.of('1 millis')
        1.sec == Duration.of('1 sec')
        2.min  == Duration.of('2 min')
        3.hour == Duration.of('3 hours')
        4.days == Duration.of('4 days')
    }

    def 'should multiply number with a duration unit' () {

        expect:
        2.sec * 5 == Duration.of('10 sec')
        5 * 3.sec == Duration.of('15 sec')

    }

    def 'should multiply number with a memory unit'() {
        expect:
        100.KB * 2 == MemoryUnit.of('200 KB')
        10 * 100.KB == MemoryUnit.of('1000 KB')
    }


    def 'should parse memory units' () {
        expect:
        1.kb  == MemoryUnit.of('1 KB')
        2.mb == MemoryUnit.of('2 MB')
        3.gb  == MemoryUnit.of('3 GB')
        4.tb == MemoryUnit.of('4 TB')
        5.pb == MemoryUnit.of('5 PB')

        3.5.GB == MemoryUnit.of('3.5 GB')
    }


}
