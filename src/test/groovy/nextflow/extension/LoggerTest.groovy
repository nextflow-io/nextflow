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

import groovy.util.logging.Slf4j
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LoggerTest extends Specification {

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @Rule
    OutputCapture capture = new OutputCapture()


    def 'should print only a debug line in a second'() {

        when:
        Dummy.log.debug1("Hello world", throttle: '1s')
        Dummy.log.debug1("Hello world", throttle: '1s')
        Dummy.log.debug1("Hello world", throttle: '1s')
        sleep (1_100)
        Dummy.log.debug1("Hello world", throttle: '1s')

        then:
        capture.toString().readLines().findAll { it.endsWith('Hello world') }.size() ==2
    }


    def 'should changed lines'() {

        when:
        Dummy.log.debug1("Hello world 1", throttle: '1s')
        Dummy.log.debug1("Hello world 2", throttle: '1s')
        Dummy.log.debug1("Hello world 1", throttle: '1s')
        Dummy.log.debug1("Hello world 2", throttle: '1s')

        then:
        capture.toString().readLines().findAll { it.contains('Hello world') }.size() ==2
        capture.toString().readLines().findAll { it.endsWith('Hello world 1') }.size() ==1
        capture.toString().readLines().findAll { it.endsWith('Hello world 2') }.size() ==1
    }


    @Slf4j
    static class Dummy {
    }


}
