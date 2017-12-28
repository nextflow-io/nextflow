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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class VersionNumberTest extends Specification {

    @Unroll
    def 'should parse version number: #str' () {

        given:
        def version = new VersionNumber(num)
        expect:
        version.major == major
        version.minor == minor
        version.patch == patch
        version.toString() == str

        where:
        num         | major     | minor     | patch | str
        null        | '0'       | null      | null  | '0'
        '0.1'       | '0'       | '1'       | null  | '0.1'
        '0.1.'      | '0'       | '1'       | null  | '0.1'
        '0.1.2'     | '0'       | '1'       | '2'   | '0.1.2'
        '3.2.1.6.4' | '3'       | '2'       | '1'   | '3.2.1.6.4'
        '1.2.RC1'   | '1'       | '2'       | 'RC1' | '1.2.RC1'

    }

    def 'should access version components' () {
        expect:
        new VersionNumber('3.5.CR2')[0] == '3'
        new VersionNumber('3.5.CR2')[1] == '5'
        new VersionNumber('3.5.CR2')[2] == 'CR2'
    }


    def 'should compare version' () {

        expect:
        new VersionNumber('1') > new VersionNumber('0')
        new VersionNumber('10') > new VersionNumber('2.1')
        new VersionNumber('1.1') < new VersionNumber('2.1')
        new VersionNumber('1.1') == new VersionNumber('1.1')
        new VersionNumber('1.1') == new VersionNumber('1.1.0')
        new VersionNumber('2.3.127') > new VersionNumber('2.3.3')
        new VersionNumber('2.3.1a') < new VersionNumber('2.3.1b')
        new VersionNumber('2.0') == new VersionNumber('2.0.0')
        new VersionNumber('2.0') < new VersionNumber('2.0.1')
        new VersionNumber('3') >  new VersionNumber('2.0.1')
        new VersionNumber('1.1') > null
        new VersionNumber('1.1') < '1.2'
    }

    @Unroll
    def 'should validate version check for version #version' () {

        expect:
        new VersionNumber(version) .matches(condition) == expected

        where:
        version     | condition     | expected
        '1.2'       | '> 1.0'       | true
        '1.2'       | '>= 1.0'      | true
        '1.2'       | '!= 1.0'      | true
        '1.2'       | '< 2.1'       | true
        '1.2'       | '<= 2.1'      | true
        '1.2'       | '<> 2.1'      | true
        '1.2'       | '1.2'         | true

        '1.2'       | '= 1.0'       | false
        '1.2'       | '< 1.0'       | false

    }

    @Unroll
    def 'should validate plus syntax for version: #version' () {

        expect:
        new VersionNumber(version) .matches(condition) == expected

        where:
        version     | condition     | expected
        '1.0'       | '1.0'         | true
        '1.0.1'     | '1.0'         | false
        '1.1'       | '1.0'         | false
        '1.2'       | '1.2.0'       | true
        '1.2'       | '1.2.1'       | false

        '1.2'       | '1.0+'        | true
        '1.2'       | '1.2+'        | true
        '1.2'       | '1.3+'        | false
        '1.2'       | '1.2.+'       | true
        '1.3'       | '1.2.+'       | false
        '1.2'       | '1.2.1+'      | false
        '1.2.1'     | '1.2.2+'      | false
        '1.2'       | '1.3+'        | false
        '1.2'       | '1+'          | true
        '1.2'       | '2+'          | false

        '2.3.1b'    | '2.3.+'       | true
        '2.4'       | '2.3.+'       | false
        '2.4'       | '2.3+'        | true
        '2.3.1b'    | '2.3.1a+'     | true
        '2.3.1b'    | '2.3.1c+'     | false
        '2.3.1b'    | '2.3.1+'      | true
        '2.3.1b'    | '2.3.1c+'     | false

        '1.2'       | '1+'          | true
        '1.2'       | '1.+'         | true
        '1.2'       | '2+'          | false
    }

}
