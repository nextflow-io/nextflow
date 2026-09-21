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

import nextflow.SysEnv
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import spock.lang.Requires
import spock.lang.Specification

/**
 * End-to-end test exercising M5 (tools + structured output together) against the real
 * OpenAI integration: the model is asked to uppercase a word, which it can only do by
 * calling the advertised {@code uppercase} tool; the tool loop runs schema-free, then a
 * final structuring turn re-encodes the free-text answer into JSON matching the declared
 * {@code {result:string}} schema. Skipped automatically when OPENAI_API_KEY is not set.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({ System.getenv('OPENAI_API_KEY') })
class AgentToolStructuredEndToEndTest extends Specification {

    private static final Map UPPERCASE_INPUT_SCHEMA = [
        type: 'object',
        properties: [text: [type: 'string']],
        required: ['text'],
        additionalProperties: false,
    ]

    private static final Map RESULT_OUTPUT_SCHEMA = [
        type: 'object',
        properties: [result: [type: 'string']],
        required: ['result'],
        additionalProperties: false,
    ]

    def 'should drive the tool loop then structure the final answer end-to-end against OpenAI'() {
        given: 'a dispatcher that uppercases the text arg and records that it ran'
        boolean dispatched = false
        ToolDispatcher dispatch = { String name, String argsJson ->
            dispatched = true
            def args = new JsonSlurper().parseText(argsJson) as Map
            return JsonOutput.toJson([result: (args.text as String).toUpperCase()])
        } as ToolDispatcher

        and: 'a request advertising the uppercase tool AND a structured {result:string} output'
        def descriptor = new ToolDescriptor(
            'uppercase', 'Uppercase the given text', UPPERCASE_INPUT_SCHEMA, [:])
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            instruction: 'To uppercase text you MUST call the uppercase tool, then reply with the tool\'s result.',
            prompt: 'uppercase the word hello',
            maxIterations: 5, tools: [], outputSchema: RESULT_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, agentName: 'uppercaser',
            // the runner never reads the environment: core resolves the credential and
            // carries it on the request, so a direct-to-runner test must supply it -- through
            // SysEnv, the same seam the production ladder uses, so a SysEnv.push here would
            // actually take effect
            apiKey: SysEnv.get('OPENAI_API_KEY'))

        when:
        def answer = new LangChainAgentRunner().run(req)

        then: 'the tool was actually called'
        dispatched

        and: 'the final answer is schema-valid JSON whose result is the uppercased word'
        def parsed = new JsonSlurper().parseText(answer) as Map
        parsed.containsKey('result')
        parsed.result.toString().toUpperCase().contains('HELLO')
    }
}
