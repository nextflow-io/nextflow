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

import dev.langchain4j.agent.tool.ToolExecutionRequest
import dev.langchain4j.agent.tool.ToolSpecification
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.invocation.InvocationContext
import dev.langchain4j.service.tool.ToolExecutionResult
import dev.langchain4j.service.tool.ToolExecutor
import dev.langchain4j.service.tool.ToolProvider
import dev.langchain4j.service.tool.ToolProviderRequest
import dev.langchain4j.service.tool.ToolProviderResult
import spock.lang.Specification

/**
 * The tracing wrapper around a skills tool provider must (a) report each skill-tool
 * invocation to the {@link AgentTrace} like a module tool, and (b) delegate
 * {@code executeWithContext} (langchain4j-skills' real execution path) so activation
 * is not broken.
 */
class SkillTraceProviderTest extends Specification {

    def 'should report skill-tool calls to the trace while delegating executeWithContext'() {
        given: 'a delegate provider exposing one skill tool whose executeWithContext returns a result'
        def spec = ToolSpecification.builder().name('activate_skill').build()
        boolean delegated = false
        def innerExec = new ToolExecutor() {
            String execute(ToolExecutionRequest req, Object memoryId) { 'via-execute' }
            ToolExecutionResult executeWithContext(ToolExecutionRequest req, InvocationContext ctx) {
                delegated = true
                return ToolExecutionResult.builder().resultText('SKILL INSTRUCTIONS').build()
            }
        }
        def delegate = new ToolProvider() {
            ToolProviderResult provideTools(ToolProviderRequest req) {
                new ToolProviderResult([(spec): innerExec] as Map<ToolSpecification,ToolExecutor>)
            }
        }
        def trace = Mock(AgentTrace)

        when: 'the wrapped provider is asked for tools and the wrapped executor runs'
        def wrapped = LangChainAgentRunner.tracingSkillProvider(delegate, trace)
        def exec = wrapped.provideTools(new ToolProviderRequest(null, UserMessage.from('hi'))).tools()[spec]
        def ter = ToolExecutionRequest.builder().id('1').name('activate_skill').arguments('{"skill_name":"greet"}').build()
        def result = exec.executeWithContext(ter, null)

        then: 'the real executor was delegated to and its result returned unchanged'
        delegated
        result.resultText() == 'SKILL INSTRUCTIONS'

        and: 'the invocation was itemized on the trace, like a module tool'
        1 * trace.tool('activate_skill', '{"skill_name":"greet"}', 'SKILL INSTRUCTIONS')
    }
}
