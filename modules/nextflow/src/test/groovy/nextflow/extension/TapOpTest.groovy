/*
 * Copyright 2013-2024, Seqera Labs
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
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class TapOpTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def 'should `tap` item to a new channel' () {
        when:
        def result = Channel.of( 4,7,9 ) .tap { first }.map { it+1 }
        then:
        session.binding.first.unwrap() == 4
        session.binding.first.unwrap() == 7
        session.binding.first.unwrap() == 9
        session.binding.first.unwrap() == Channel.STOP

        result.unwrap() == 5
        result.unwrap() == 8
        result.unwrap() == 10
        result.unwrap() == Channel.STOP

        !session.dag.isEmpty()
    }

    def 'should `tap` item to more than one channel' () {
        when:
        def result = Channel.of( 4,7,9 ) .tap { foo; bar }.map { it+1 }
        then:
        session.binding.foo.unwrap() == 4
        session.binding.foo.unwrap() == 7
        session.binding.foo.unwrap() == 9
        session.binding.foo.unwrap() == Channel.STOP
        session.binding.bar.unwrap() == 4
        session.binding.bar.unwrap() == 7
        session.binding.bar.unwrap() == 9
        session.binding.bar.unwrap() == Channel.STOP

        result.unwrap() == 5
        result.unwrap() == 8
        result.unwrap() == 10
        result.unwrap() == Channel.STOP

        !session.dag.isEmpty()
    }

    def 'should `tap` item with data value' () {
        when:
        def result = Channel.value(1 ) .tap { first }.map { it+1 }
        then:
        session.binding.first.unwrap() == 1
        session.binding.first.unwrap() == 1
        and:
        result.unwrap() == 2
        result.unwrap() == 2

        !session.dag.isEmpty()
    }

}
