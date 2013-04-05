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

class FileSizeTest extends Specification {

    def 'test toString'() {

        expect:
        new FileSize(1024).toString() == '1 KB'
        new FileSize(2 * 1024 * 1024).toString() == '2 MB'
        new FileSize(3 * 1024L * 1024L * 1024L).toString() == '3 GB'
        new FileSize(4 * 1024L * 1024L * 1024L * 1024L).toString() == '4 TB'
        new FileSize(5 * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '5 PB'
        new FileSize(6 * 1024L * 1024L * 1024L * 1024L * 1024L * 1024L).toString() == '6 EB'

    }


    def 'test create' () {

        expect:
        new FileSize('1').toBytes() == 1
        new FileSize('2B').toBytes() == 2
        new FileSize('3 B').toBytes() == 3

        new FileSize('1KB').toBytes() == 1024
        new FileSize('2 KB').toBytes() == 2 * 1024
        new FileSize('3.5 KB').toBytes() == 3.5 * 1024

        new FileSize('1MB').toBytes() == 1024 * 1024
        new FileSize('2 MB').toBytes() == 2 * 1024 * 1024
        new FileSize('3.5 MB').toBytes() == 3.5 * 1024 * 1024

        new FileSize('1GB').toBytes() == 1024 * 1024 * 1024L
        new FileSize('2 GB').toBytes() == 2 * 1024 * 1024 * 1024L
        new FileSize('3.5 GB').toBytes() == 3.5 * 1024 * 1024 * 1024L

        new FileSize('1TB').toBytes() == 1024 * 1024 * 1024L* 1024L
        new FileSize('2 TB').toBytes() == 2 * 1024 * 1024 * 1024L* 1024L
        new FileSize('3.5 TB').toBytes() == 3.5 * 1024 * 1024 * 1024L* 1024L

        new FileSize('1PB').toBytes() == 1024 * 1024 * 1024L* 1024L* 1024L
        new FileSize('2 PB').toBytes() == 2 * 1024 * 1024 * 1024L* 1024L* 1024L
        new FileSize('3.5 PB').toBytes() == 3.5 * 1024 * 1024 * 1024L* 1024L* 1024L
    }

}