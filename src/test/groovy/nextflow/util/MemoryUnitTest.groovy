/*
 * Copyright (c) 2012, the authors.
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

}