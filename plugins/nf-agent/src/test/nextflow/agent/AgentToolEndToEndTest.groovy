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

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import spock.lang.Requires
import spock.lang.Specification

/**
 * End-to-end test exercising the FULL tool-call loop through the real OpenAI
 * integration: the model is asked to uppercase a word, which it can only do by
 * calling the {@code uppercase} tool we advertise; the runner dispatches the
 * call back to our {@link ToolDispatcher} and feeds the result back to the model
 * until it returns a final answer. Skipped automatically when OPENAI_API_KEY is
 * not set (CI and keyless dev environments).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({ System.getenv('OPENAI_API_KEY') })
class AgentToolEndToEndTest extends Specification {

    private static final Map UPPERCASE_INPUT_SCHEMA = [
        type: 'object',
        properties: [text: [type: 'string']],
        required: ['text'],
        additionalProperties: false,
    ]

    def 'should drive the real tool-call loop end-to-end against OpenAI'() {
        given: 'a dispatcher that uppercases the text arg and records that it ran'
        boolean dispatched = false
        ToolDispatcher dispatch = { String name, String argsJson ->
            dispatched = true
            def args = new JsonSlurper().parseText(argsJson) as Map
            def text = args.text as String
            return JsonOutput.toJson([result: text.toUpperCase()])
        } as ToolDispatcher

        and: 'a request advertising the single uppercase tool'
        def descriptor = new ToolDescriptor(
            'uppercase',
            'Uppercase the given text',
            UPPERCASE_INPUT_SCHEMA,
            [:])
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'To uppercase text you MUST call the uppercase tool, then reply with only the tool\'s result.',
            'uppercase the word hello',
            5,
            [],
            null,
            null,
            [descriptor],
            dispatch)

        when:
        def answer = new LangChainAgentRunner().run(req)

        then: 'the tool was actually called'
        dispatched

        and: 'the final answer reflects the uppercased result'
        answer.toUpperCase().contains('HELLO')
    }
}
