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

import spock.lang.Specification

class AgentRunnerRequestTest extends Specification {

    def 'should build via named args with goal as the last field'() {
        when:
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            instruction: 'sys',
            prompt: 'p',
            maxIterations: 7,
            tools: [],
            outputSchema: null,
            inputJson: '{}',
            toolSpecs: null,
            dispatch: null,
            requestTimeoutSeconds: 30,
            goal: 'reach the objective')

        then:
        req.model == 'openai/gpt-5-mini'
        req.instruction == 'sys'
        req.prompt == 'p'
        req.maxIterations == 7
        req.inputJson == '{}'
        req.requestTimeoutSeconds == 30
        req.goal == 'reach the objective'
    }

    def 'should default goal to null when omitted'() {
        when:
        def req = new AgentRunnerRequest(model: 'm', prompt: 'p')
        then:
        req.goal == null
    }
}
