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

package nextflow.extension
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowChoiceExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }


    @Timeout(1)
    def 'should choice targets specified with an open array'() {

        given:
        def source = Channel.from 'Hello world', 'Hola', 'Hello John'
        def queue1 = Channel.create()
        def queue2 = Channel.create()
        def queue3 = Channel.create()

        when:
        source.choice( queue1, queue2, queue3 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue1.val == 'Hello John'
        queue1.val == Channel.STOP
        queue2.val == 'Hola'
        queue2.val == Channel.STOP
        queue3.val == Channel.STOP

    }

    def 'should choice targets specified with a list '() {

        given:
        def source = Channel.from 'Hello world', 'Hola', 'Hello John'
        def queue1 = Channel.create()
        def queue2 = Channel.create()

        when:
        source.choice( [queue1, queue2] ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue1.val == 'Hello John'
        queue1.val == Channel.STOP
        queue2.val == 'Hola'
        queue2.val == Channel.STOP

    }


    @Timeout(1)
    def 'should choice from a dataflow value channel'() {

        given:
        def source = Channel.value('Hello world')
        def queue1 = new DataflowVariable()
        def queue2 = new DataflowVariable()

        when:
        source.choice( queue1, queue2 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue2.val == Channel.STOP

    }


}
