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

        '1.2'       | '>1.0, <=1.2' | true
        '1.2'       | '>1.0, <1.1'  | false

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

    def 'should match hashCode' () {
        expect:
        new VersionNumber('1.0.0').hashCode() == new VersionNumber('1.0.0').hashCode()
        new VersionNumber('1.0.0').hashCode() != new VersionNumber('1.0.1').hashCode()
    }
}
