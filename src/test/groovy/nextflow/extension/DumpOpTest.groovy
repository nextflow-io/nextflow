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

package nextflow.extension

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll
import test.OutputCapture
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DumpOpTest extends Specification {

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    def 'should dump channel items'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.from(1,2,3).dump()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP] 1')
        capture.toString().contains('[DUMP] 2')
        capture.toString().contains('[DUMP] 3')

    }

    def 'should dump channel items with closure'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.from(1,2,3).dump { it * it }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP] 1')
        capture.toString().contains('[DUMP] 4')
        capture.toString().contains('[DUMP] 9')

    }

    def 'should print channel items with a tag'() {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = Channel.from(1,2,3).dump(tag:'foo')
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString().contains('[DUMP: foo] 1')
        capture.toString().contains('[DUMP: foo] 2')
        capture.toString().contains('[DUMP: foo] 3')

    }

    @Unroll
    def 'should validate isEnabled when tag=#tag and names=#names' () {

        given:
        def op = new DumpOp(tag: tag, dumpNames: names.tokenize(','))

        expect:
        op.isEnabled() == expected

        where:
        tag     | names         | expected
        'foo'   | 'foo'         | true
        'foo'   | '*'           | true
        'foo'   | 'bar'         | false
        'foo'   | 'f'           | false
        'foo'   | 'f*'          | true
        'bar'   | 'f*'          | false
        'bar'   | 'f*,b*'       | true
        null    | '*'           | true

    }

}
