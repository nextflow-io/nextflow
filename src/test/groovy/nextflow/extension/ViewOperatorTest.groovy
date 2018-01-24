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
import spock.lang.Specification
import test.OutputCapture

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ViewOperatorTest extends Specification{

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @org.junit.Rule
    OutputCapture capture = new OutputCapture()


    def 'should print channel items'() {

        when:
        def result = Channel.from(1,2,3).view()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString() == '1\n2\n3\n'

    }

    def 'should print channel items applying the closure formatting rule'() {

        when:
        def result = Channel.from(1,2,3).view { "~ $it " }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        capture.toString() == '~ 1 \n~ 2 \n~ 3 \n'

    }

    def 'should print channel items without appending the newline character'() {

        when:
        def result = Channel.from(1,2,3).view(newLine:false) { " ~ $it" }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        capture.toString() == ' ~ 1 ~ 2 ~ 3'
    }

    def 'should print dataflow value' () {
        when:
        def result = Channel.value(1).view { ">> $it" }
        then:
        result.val == 1
        capture.toString() == ">> 1\n"
    }

    def 'should return stop signal'() {
        when:
        def result = Channel.value(Channel.STOP).view { ">> $it" }
        then:
        result.val == Channel.STOP
        capture.toString() == ''
    }

}
