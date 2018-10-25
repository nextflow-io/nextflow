/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
