/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.util

import java.util.concurrent.CompletableFuture
import java.util.concurrent.TimeUnit

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SimpleAgentTest extends Specification {

    def 'should update agent' () {
        given:
        def state = []
        def agent = new SimpleAgent(state)

        when:
        agent.send { state<<1 }
        agent.send { state<<2 }
        agent.send { state<<3 }
        and:
        def result = agent.getValue()

        then:
        result == [1,2,3]
        and:
        agent.getValue() == agent.getQuickValue()

        // changing the state object does not modify the
        // result because it's a clone of the original one
        when:
        state << 4
        then:
        result == [1,2,3]

    }

    def 'should get the value when invoked by the runner thread' () {
        given:
        def state = []
        def result = new CompletableFuture()
        SimpleAgent agent
        // the error handler is invoked by the agent runner thread itself, therefore
        // it cannot await for an event that only that thread is able to serve
        agent = new SimpleAgent(state).onError { result.complete(agent.getValue()) }

        when:
        agent.send { state<<1 }
        agent.send { throw new RuntimeException('Oops') }

        then:
        result.get(30, TimeUnit.SECONDS) == [1]
    }

}
