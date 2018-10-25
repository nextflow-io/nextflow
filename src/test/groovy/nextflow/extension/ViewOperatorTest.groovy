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
