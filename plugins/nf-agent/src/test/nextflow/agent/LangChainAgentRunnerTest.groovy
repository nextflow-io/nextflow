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

import dev.langchain4j.data.message.AiMessage
import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.response.ChatResponse
import spock.lang.Specification

class LangChainAgentRunnerTest extends Specification {

    def 'should send instruction + prompt and return the assistant text'() {
        given:
        List<ChatMessage> captured = null
        // langchain4j ChatModel has no single abstract method (all default), so
        // mock it by overriding the chat(List<ChatMessage>) entry point used by the runner.
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                captured = messages
                ChatResponse.builder().aiMessage(AiMessage.from('CANNED ANSWER')).build()
            }
        ] as ChatModel

        and:
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(_, _) >> model
        })
        def req = new AgentRunnerRequest('openai/gpt-5-mini', 'You are helpful.', 'Question: hi', 5, [])

        when:
        def answer = runner.run(req)

        then:
        answer == 'CANNED ANSWER'
        and:
        captured.size() == 2
        captured[0] instanceof SystemMessage
        (captured[0] as SystemMessage).text() == 'You are helpful.'
        captured[1] instanceof UserMessage
        (captured[1] as UserMessage).singleText() == 'Question: hi'
    }

    def 'should omit the system message when no instruction is given'() {
        given:
        List<ChatMessage> captured = null
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                captured = messages
                ChatResponse.builder().aiMessage(AiMessage.from('OK')).build()
            }
        ] as ChatModel

        and:
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(_, _) >> model
        })
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'just a prompt', 5, [])

        when:
        def answer = runner.run(req)

        then:
        answer == 'OK'
        captured.size() == 1
        captured[0] instanceof UserMessage
    }

    def 'should fail when the model is missing'() {
        given:
        def runner = new LangChainAgentRunner()
        def req = new AgentRunnerRequest(null, 'inst', 'prompt', 5, [])

        when:
        runner.run(req)

        then:
        thrown(IllegalArgumentException)
    }
}
