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

import nextflow.Session
import spock.lang.Specification

class AgentFactoryTest extends Specification {

    def 'should build an AgentDef from a body closure'() {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def factory = new AgentFactory(script, session)

        when:
        def agent = factory.newAgent('eval_agent', { -> })

        then:
        agent instanceof AgentDef
        agent.name == 'eval_agent'
    }
}
