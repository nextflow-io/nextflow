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
import dev.langchain4j.model.chat.listener.ChatModelListener
import dev.langchain4j.model.chat.listener.ChatModelResponseContext
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Renders a readable trace of an agent run that simulates the agent "thinking": each
 * model round-trip is a turn, the model's free text is shown as its reasoning, and each
 * tool invocation is reported as it runs.
 *
 * <p>Two levels: the human-readable narrative — turns, model reasoning, which tool was
 * invoked, the final answer — is logged at INFO; the low-level tool inputs and outputs
 * are logged at DEBUG.
 *
 * <p>Turn boundaries and reasoning are captured as a langchain4j
 * {@link ChatModelListener} ({@link #onResponse} fires once per model round-trip); tool
 * invocations are reported via {@link #tool} from the runner's tool executor. Everything
 * runs on the single agent operator thread (AiServices dispatches sequentially), so the
 * turn counter needs no synchronization.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AgentTrace implements ChatModelListener {

    private final String label

    private int turn

    AgentTrace(String agentName) {
        this.label = agentName ? "[agent ${agentName}]" : '[agent]'
    }

    /** Emit the run header: model and the tools the agent may call. */
    void begin(String model, Collection<String> toolNames) {
        emit("model ${model}" + (toolNames ? " · tools: ${toolNames.join(', ')}" : ' · no tools'))
    }

    /** Emit the run footer with the final answer. */
    void end(String finalAnswer) {
        paragraph('final answer: ', finalAnswer)
    }

    @Override
    void onResponse(ChatModelResponseContext ctx) {
        final AiMessage msg = ctx.chatResponse().aiMessage()
        turn += 1
        final boolean isFinal = !msg.hasToolExecutionRequests()
        emit(isFinal ? "── turn ${turn} ── (final)" : "── turn ${turn} ──")
        // the model's reasoning: its dedicated thinking channel when present, else its free-text
        // content on a tool-decision turn (on the final turn the content is the answer, shown by
        // end()). NOTE: most OpenAI tool-call turns carry no content/thinking, so a reasoning line
        // appears only when the model actually narrates (or returnThinking is supported by the model).
        final String reasoning = msg.thinking()?.trim() ?: (isFinal ? null : msg.text()?.trim())
        if( reasoning )
            paragraph('reasoning: ', reasoning)
        // tool inputs/outputs are emitted by #tool as each call runs
    }

    /**
     * Report a single tool invocation: the decision — the tool name plus a short, readable digest
     * of its arguments (paths shown as file names, nested maps inlined) — is the human-readable
     * narrative logged at INFO; the full raw JSON arguments and result are logged at DEBUG.
     */
    void tool(String name, String argsJson, String resultJson) {
        emit("  → ${name}(${summarizeArgs(argsJson)})")
        emitDebug("  ${name} input: ${argsJson}")
        emitDebug("  ${name} output: ${resultJson}")
    }

    private void emit(String text) {
        log.info("${label} ${text}")
    }

    private void emitDebug(String text) {
        log.debug("${label} ${text}")
    }

    /**
     * Emit possibly multi-line text under a marker: the first line carries the
     * marker, continuation lines are indented so the block reads as one thought.
     */
    private void paragraph(String marker, String text) {
        final lines = (text ?: '').readLines()
        if( !lines ) {
            emit(marker.trim())
            return
        }
        emit("  ${marker}${lines[0]}")
        final pad = ' ' * (marker.length() + 2)
        for( int i = 1; i < lines.size(); i++ )
            emit("${pad}${lines[i]}")
    }

    /** Max length of the rendered argument digest on the INFO line. */
    private static final int MAX_ARGS_CHARS = 100

    /**
     * A short, human-readable digest of a tool's JSON arguments for the INFO line: top-level
     * {@code key=value} pairs, with file-path values shown as their file name, nested maps inlined
     * as {@code {k:v}}, and the whole thing clipped to {@link #MAX_ARGS_CHARS}. Falls back to the
     * flattened raw string when the arguments are not a JSON object.
     */
    private static String summarizeArgs(String argsJson) {
        final parsed = parseJson(argsJson)
        if( !(parsed instanceof Map) )
            return clip(flatten(argsJson), MAX_ARGS_CHARS)
        final entries = ((Map) parsed).collect { k, v -> "${k}=${renderValue(v)}".toString() }
        return clip(entries.join(', '), MAX_ARGS_CHARS)
    }

    private static String renderValue(Object v) {
        if( v == null )
            return 'null'
        if( v instanceof Map ) {
            final inner = ((Map) v).collect { k, x -> "${k}:${scalar(x)}".toString() }.join(', ')
            return "{${inner}}"
        }
        if( v instanceof List ) {
            final list = (List) v
            return list.isEmpty() ? '[]' : "[${scalar(list[0])}${list.size() > 1 ? ', …' : ''}]"
        }
        return scalar(v)
    }

    /** Render a scalar: a path-like value becomes its file name; anything long is clipped. */
    private static String scalar(Object v) {
        if( v == null )
            return 'null'
        final s = v.toString()
        final slash = s.lastIndexOf('/')
        return slash >= 0 ? s.substring(slash + 1) : clip(s, 40)
    }

    private static String clip(String s, int max) {
        if( s == null )
            return ''
        return s.length() > max ? s.substring(0, max) + '…' : s
    }

    private static String flatten(String s) {
        return s != null ? s.replaceAll(/\s+/, ' ').trim() : ''
    }

    private static Object parseJson(String json) {
        if( !json?.trim() )
            return null
        try {
            return new JsonSlurper().parseText(json)
        }
        catch( Exception e ) {
            log.trace("Agent trace: tool arguments are not valid JSON, showing raw: ${e.message}")
            return null
        }
    }
}
