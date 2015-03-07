/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import test.OutputCapture
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowExtensionsTest2 extends Specification{

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    def testPrint() {

        when:
        Channel.from(1,2,3).print()
        sleep 50
        then:
        capture.toString() == '123'

    }

    def testPrintWithClosure() {

        when:
        Channel.from(1,2,3).print { "> $it " }
        sleep 50
        then:
        capture.toString() == '> 1 > 2 > 3 '

    }

    def testPrintln() {

        when:
        Channel.from(1,2,3).println()
        sleep 50
        then:
        capture.toString() == '1\n2\n3\n'

    }

    def testPrintlnWithClosure() {

        when:
        Channel.from(1,2,3).map{ "> $it" }.println()
        sleep 50
        then:
        capture.toString() == '> 1\n> 2\n> 3\n'

    }


    def testView() {

        when:
        def result = Channel.from(1,2,3).view()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        capture.toString() == '1\n2\n3\n'

    }

    def testViewWithClosure() {

        when:
        def result = Channel.from(1,2,3).view { "~ $it " }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        capture.toString() == '~ 1 \n~ 2 \n~ 3 \n'

    }

    def testViewNoNewLine() {

        when:
        def result = Channel.from(1,2,3).view(newLine:false) { " ~ $it" }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

        capture.toString() == ' ~ 1 ~ 2 ~ 3'
    }

}
