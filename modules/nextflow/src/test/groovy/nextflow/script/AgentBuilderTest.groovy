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

class AgentBuilderTest extends Specification {

    def 'should capture directives, inputs, outputs and prompt then build an AgentDef'() {
        given:
        def builder = new AgentBuilder(null, 'eval_agent')

        when:
        builder.model('openai/gpt-5-mini')
        builder.instruction('You are helpful.')
        builder.tools()
        builder.maxIterations(20)
        builder._input_('question', String)
        builder._output_('plan', String)
        def prompt = new PromptDef({ "Question: x" }, 'Question: ${question}')
        def agent = builder.withPrompt(prompt).build()

        then:
        agent instanceof AgentDef
        agent.name == 'eval_agent'
        agent.model == 'openai/gpt-5-mini'
        agent.instruction == 'You are helpful.'
        agent.maxIterations == 20
        agent.tools == []
        agent.inputs*.name == ['question']
        agent.outputs*.name == ['plan']
        agent.prompt.source == 'Question: ${question}'
    }

    def 'should capture the skills directive (single and list)'() {
        given:
        def builder = new AgentBuilder(null, 'a')
        def b2 = new AgentBuilder(null, 'b')

        when:
        builder.skills('greet')
        builder._input_('q', String); builder._output_('a', String)
        def a1 = builder.withPrompt(new PromptDef({ 'x' }, 'x')).build()
        b2.skills('greet', 'github.com/org/repo')
        b2._input_('q', String); b2._output_('a', String)
        def a2 = b2.withPrompt(new PromptDef({ 'x' }, 'x')).build()

        then:
        a1.skills == ['greet']
        a2.skills == ['greet', 'github.com/org/repo']
    }

    def 'should reject an unknown directive'() {
        given:
        def builder = new AgentBuilder(null, 'a')

        when:
        builder.bogusDirective('x')

        then:
        thrown(Exception)
    }

    def 'should fail to build without a prompt'() {
        given:
        def builder = new AgentBuilder(null, 'a')

        when:
        builder.build()

        then:
        thrown(IllegalStateException)
    }
}
