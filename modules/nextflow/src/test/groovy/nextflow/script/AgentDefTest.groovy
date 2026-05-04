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
package nextflow.script

import spock.lang.Specification

class AgentDefTest extends Specification {

    def 'should construct an AgentDef with name and body'() {
        given:
        def script = Mock(BaseScript)
        def body = { -> /* placeholder */ } as Closure

        when:
        def agent = new AgentDef(script, body, 'eval_agent')

        then:
        agent.name == 'eval_agent'
        agent.simpleName == 'eval_agent'
        agent.type == 'agent'
    }

    def 'should throw on run() since execution is not yet implemented'() {
        given:
        def script = Mock(BaseScript)
        def body = { -> } as Closure
        def agent = new AgentDef(script, body, 'foo')

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(UnsupportedOperationException)
        e.message.contains('agent execution not yet implemented')
    }

    def 'should clone with a new name'() {
        given:
        def script = Mock(BaseScript)
        def agent = new AgentDef(script, { -> }, 'foo')

        when:
        def renamed = agent.cloneWithName('bar')

        then:
        renamed instanceof AgentDef
        renamed.name == 'bar'
        agent.name == 'foo' // original untouched
    }
}
