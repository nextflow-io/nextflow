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

import groovy.json.JsonSlurper
import spock.lang.Requires
import spock.lang.Specification

/**
 * End-to-end test exercising the real OpenAI integration through the
 * langchain4j runner, using structured (JSON-schema) output. Skipped
 * automatically when OPENAI_API_KEY is not set (CI and keyless dev
 * environments).
 */
@Requires({ System.getenv('OPENAI_API_KEY') })
class AgentEndToEndTest extends Specification {

    def 'should get a real structured answer from the model'() {
        given:
        def runner = new LangChainAgentRunner()
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'You are a terse geography assistant.',
            'What is the capital of France?',
            5,
            [],
            [
                type: 'object',
                properties: [
                    capital: [type: 'string'],
                    country: [type: 'string'],
                ],
                required: ['capital', 'country'],
                additionalProperties: false,
            ],
            null)

        when:
        def answer = runner.run(req)

        then:
        answer != null
        !answer.trim().isEmpty()
        and:
        // the model returns structured JSON matching the requested schema
        def parsed = new JsonSlurper().parseText(answer) as Map
        parsed.capital.toString().toLowerCase().contains('paris')
    }
}
