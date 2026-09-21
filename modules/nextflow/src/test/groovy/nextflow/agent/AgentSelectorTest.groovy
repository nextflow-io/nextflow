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
package nextflow.agent

import nextflow.Global
import nextflow.Session
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.AgentDef
import nextflow.script.BaseScript
import nextflow.script.PromptDef
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockSession

/**
 * Aliasing side of the agent config-selector story (module-include design §6.1/§6.2, as amended by
 * its §16 -- selectors live in the independent {@code agent} scope, not {@code process}). The
 * premise of a module is that its CONSUMER cannot edit the module file: placing an included agent
 * from the config must therefore work for the DECLARED name as well as the alias -- which is why
 * {@link AgentDef} carries a stable {@code baseName} like {@code ProcessDef} does.
 *
 * <p>The selector resolution itself (the precedence ladder, the `agent` vs `process` scope split)
 * is covered by {@code AgentConfigSelectorTest}; this spec covers what {@code cloneWithName}
 * preserves and what it registers.
 *
 * <p>Extends {@code Dsl2Spec} because {@code cloneWithName} mutates the static
 * {@code ScriptMeta.resolvedAgentNames} and {@code Dsl2Spec.setup()} resets it.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
class AgentSelectorTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        AgentCallInfo.clear()
    }

    private Session newSession(Map config = null) {
        final session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()
        Global.session = session
        return session
    }

    /** An agent DECLARED as `qa`, as an agent module would declare it. */
    private AgentDef newAgent(Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    def 'cloneWithName changes the name and the simple name but keeps the base name'() {
        given:
        final agent = newAgent()

        when:
        final alias = (AgentDef) agent.cloneWithName('hi')

        then:
        alias.getName() == 'hi'
        alias.getSimpleName() == 'hi'
        alias.getBaseName() == 'qa'
        and: 'the template is untouched'
        agent.getName() == 'qa'
        agent.getBaseName() == 'qa'

        when: 'the clone name carries a workflow scope prefix'
        final scoped = (AgentDef) agent.cloneWithName('wf:a')

        then:
        scoped.getName() == 'wf:a'
        scoped.getSimpleName() == 'a'
        scoped.getBaseName() == 'qa'
    }

    // -- R2 (cosmetic): the alias and the workflow-scoped name must reach the AGENT selector
    //    registry so Session#checkConfig does not warn "no agent matching config selector" for a
    //    selector that actually applied.
    def 'the alias reaches the config-selector registry'() {
        given:
        final session = newSession([agent: ['withName:hi': [cpus: 2]]])

        when:
        newAgent().cloneWithName('hi')

        then:
        ScriptMeta.allAgentNames().contains('hi')
        and: 'so the selector is not reported as unmatched'
        session.validateConfig0([], ScriptMeta.allAgentNames()) == []
    }
}
