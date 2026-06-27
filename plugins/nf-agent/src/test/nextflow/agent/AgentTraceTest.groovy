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

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import dev.langchain4j.agent.tool.ToolExecutionRequest
import dev.langchain4j.data.message.AiMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.ModelProvider
import dev.langchain4j.model.chat.listener.ChatModelResponseContext
import dev.langchain4j.model.chat.request.ChatRequest
import dev.langchain4j.model.chat.response.ChatResponse
import org.slf4j.LoggerFactory
import spock.lang.Specification

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentTraceTest extends Specification {

    ListAppender<ILoggingEvent> appender
    Logger logger

    def setup() {
        logger = (Logger) LoggerFactory.getLogger(AgentTrace)
        logger.setLevel(Level.TRACE)
        appender = new ListAppender<ILoggingEvent>()
        appender.start()
        logger.addAppender(appender)
    }

    def cleanup() {
        logger.detachAppender(appender)
    }

    private List<String> lines() {
        return appender.list*.formattedMessage
    }

    private static ChatModelResponseContext respCtx(AiMessage msg) {
        final resp = ChatResponse.builder().aiMessage(msg).build()
        final req = ChatRequest.builder().messages(UserMessage.from('x')).build()
        return new ChatModelResponseContext(resp, req, ModelProvider.OTHER, [:])
    }

    def 'labels every line with the agent name'() {
        when:
        new AgentTrace('triage').begin('openai/gpt-5', ['FASTQC'])
        then:
        lines().every { it.startsWith('[agent triage] ') }
        lines().join('\n').contains('model openai/gpt-5')
        lines().join('\n').contains('tools: FASTQC')
    }

    def 'numbers turns and shows the model reasoning'() {
        given:
        def trace = new AgentTrace('triage')
        def decide = AiMessage.from('I will run FastQC first', [ToolExecutionRequest.builder().name('FASTQC').arguments('{}').build()])
        def finalMsg = AiMessage.from('All samples passed QC')

        when:
        trace.onResponse(respCtx(decide))   // turn 1: decides a tool
        trace.onResponse(respCtx(finalMsg)) // turn 2: final answer
        def out = lines().join('\n')

        then:
        out.contains('── turn 1 ──')
        out.contains('reasoning: I will run FastQC first')
        out.contains('── turn 2 ── (final)')
    }

    def 'logs the tool decision digest at INFO and full payloads at DEBUG'() {
        when:
        new AgentTrace('triage').tool('FASTQC', '{"reads":"/data/s1.fq.gz"}', '{"html":"/work/ab/s1.html"}')
        def infos = appender.list.findAll { it.level == Level.INFO }*.formattedMessage
        def debugs = appender.list.findAll { it.level == Level.DEBUG }*.formattedMessage

        then: 'INFO shows the decision digest (path as file name), not the raw payloads'
        infos.any { it.contains('→ FASTQC(reads=s1.fq.gz)') }
        !infos.join('\n').contains('/data/s1.fq.gz')
        !infos.join('\n').contains('/work/ab/s1.html')

        and: 'the full raw inputs and outputs are low-level DEBUG'
        debugs.any { it.contains('FASTQC input:') && it.contains('/data/s1.fq.gz') }
        debugs.any { it.contains('FASTQC output:') && it.contains('/work/ab/s1.html') }
    }

    def 'renders a readable argument digest with nested maps'() {
        when:
        new AgentTrace('a').tool('SKESA', '{"reads":"/work/ab/isolate_001.fastq.gz","meta":{"id":"isolate_001"}}', '{}')
        def infos = appender.list.findAll { it.level == Level.INFO }*.formattedMessage
        then:
        infos.any { it.contains('→ SKESA(reads=isolate_001.fastq.gz, meta={id:isolate_001})') }
    }

    def 'renders a multi-line final answer with the answer text'() {
        when:
        new AgentTrace('triage').end('line one\nline two')
        def out = lines().join('\n')
        then:
        out.contains('final answer: line one')
        out.contains('line two')
    }
}
