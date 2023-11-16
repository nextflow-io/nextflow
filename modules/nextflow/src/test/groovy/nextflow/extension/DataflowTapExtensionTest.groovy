/*
 * Copyright 2013-2023, Seqera Labs
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
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowTapExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }


    def 'should `tap` item to a new channel' () {

        when:
        def result = Channel.of( 4,7,9 ) .tap { first }.map { it+1 }
        then:
        session.binding.first.val == 4
        session.binding.first.val == 7
        session.binding.first.val == 9
        session.binding.first.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` item to more than one channel' () {

        when:
        def result = Channel.of( 4,7,9 ) .tap { foo; bar }.map { it+1 }
        then:
        session.binding.foo.val == 4
        session.binding.foo.val == 7
        session.binding.foo.val == 9
        session.binding.foo.val == Channel.STOP
        session.binding.bar.val == 4
        session.binding.bar.val == 7
        session.binding.bar.val == 9
        session.binding.bar.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

}
