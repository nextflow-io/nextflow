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
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowIntoExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }


    def 'should create two target channels'() {
        given:
        def result = Channel.from(1,2,3,4)

        when:
        def (ch1, ch2) = result.into(2)

        then:
        ch1 instanceof DataflowQueue
        ch2 instanceof DataflowQueue

        ch1.val == 1
        ch1.val == 2
        ch1.val == 3
        ch1.val == 4
        ch1.val == Channel.STOP

        ch2.val == 1
        ch2.val == 2
        ch2.val == 3
        ch2.val == 4
        ch2.val == Channel.STOP

    }

    def 'should forward items into target channels'() {
        given:
        def result = Channel.from('a','b',[1,2])
        def ch1 = Channel.create()
        def ch2 = Channel.create()
        def ch3 = Channel.create()

        when:
        result.into(ch1, ch2, ch3)

        then:
        ch1.val == 'a'
        ch1.val == 'b'
        ch1.val == [1,2]
        ch1.val == Channel.STOP

        ch2.val == 'a'
        ch2.val == 'b'
        ch2.val == [1,2]
        ch2.val == Channel.STOP

        ch3.val == 'a'
        ch3.val == 'b'
        ch3.val == [1,2]
        ch3.val == Channel.STOP

    }

    @Timeout(1)
    def 'should forward dataflow value into a new channel'() {

        when:
        def result = Channel.create()
        Channel.value('Hello').into(result)
        then:
        result.val == 'Hello'
        result.val == Channel.STOP

    }

    def 'should create new dataflow variables and forward item to them'  () {

        when:
        Channel.from(10,2,30).into { alpha; gamma }

        then:
        session.binding.alpha.val == 10
        session.binding.alpha.val == 2
        session.binding.alpha.val == 30
        session.binding.alpha.val == Channel.STOP

        session.binding.gamma.val == 10
        session.binding.gamma.val == 2
        session.binding.gamma.val == 30
        session.binding.gamma.val == Channel.STOP

    }


}
