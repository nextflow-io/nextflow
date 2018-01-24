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
import org.junit.Rule
import test.OutputCapture
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PrintOperatorTest extends Specification{

    def setupSpec() {
        new Session()
    }

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should print channel item to stdout'() {

        when:
        Channel.from(1,2,3).print()
        sleep 50
        then:
        capture.toString() == '123'

    }

    def 'should print item applying closure formatting rule'() {

        when:
        Channel.from(1,2,3).print { "> $it " }
        sleep 50
        then:
        capture.toString() == '> 1 > 2 > 3 '

    }

    def 'should print item appending newline character'() {

        when:
        Channel.from(1,2,3).println()
        sleep 50
        then:
        capture.toString() == '1\n2\n3\n'

    }

    def 'should print items applying closure formatting and appending newline'() {

        when:
        Channel.from(1,2,3).map{ "> $it" }.println()
        sleep 50
        then:
        capture.toString() == '> 1\n> 2\n> 3\n'

    }



}
