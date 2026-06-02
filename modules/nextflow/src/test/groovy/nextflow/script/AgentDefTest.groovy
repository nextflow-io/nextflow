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

import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification

class AgentDefTest extends Specification {

    private AgentDef makeAgent(BaseScript script, String name) {
        def prompt = new PromptDef({ -> 'hello' }, 'hello')
        return new AgentDef(script, name, [:], [], [], prompt)
    }

    def 'should construct an AgentDef with name and content'() {
        given:
        def script = Mock(BaseScript)

        when:
        def agent = makeAgent(script, 'eval_agent')

        then:
        agent.name == 'eval_agent'
        agent.simpleName == 'eval_agent'
        agent.type == 'agent'
    }

    def 'should fail on run() when the agent does not declare exactly one input'() {
        given:
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo') // no inputs declared

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare exactly one input')
    }

    def 'should clone with a new name'() {
        given:
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo')

        when:
        def renamed = agent.cloneWithName('bar')

        then:
        renamed instanceof AgentDef
        renamed.name == 'bar'
        agent.name == 'foo' // original untouched
    }
}
