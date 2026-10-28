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

import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type
import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic

import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHelper
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.util.TypeHelper

/**
 * Output schema and decoding strategy for an agent invocation: the mode the answer comes back in,
 * the JSON schema the model is given for it (null when there is none), and the two ways that
 * answer is turned back into the declared outputs -- {@link #decode} for a canonical task's
 * terminal frame, {@link #bind} for an in-JVM runner's return value.
 *
 * <p>Owning both keeps the mode and its consumers together: the plan is decided once while
 * lowering an agent ({@code AgentDef.resolveOutputPlan}) and then travels to the task body, where
 * it is the only thing that knows how to read the answer.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentOutputPlan {

    final AgentOutputMode mode
    final Map schema

    AgentOutputPlan(AgentOutputMode mode, Map schema) {
        this.mode = mode
        this.schema = schema
    }

    /**
     * Decode a canonical terminal frame, including a scalar final_answer wrapper.
     *
     * <p>Split where the error messages divide: reading the ANSWER out of the task's last stdout
     * frame fails with a frame-level message, and interpreting that answer per {@link #mode} fails
     * with an output-level one.
     */
    Object decode(Object stdout, String outputName, Type outputType, Path workDir = null) {
        final String answer = terminalAnswer(stdout)
        if( mode == AgentOutputMode.TEXT )
            return TypeHelper.asType(answer, outputType)
        final Object value = new JsonSlurper().parseText(stripFences(answer))
        if( mode == AgentOutputMode.SCALAR_CONTRACT ) {
            if( !(value instanceof Map) || !((Map)value).containsKey(outputName) )
                throw new ScriptRuntimeException('Canonical agent scalar output must be a JSON object containing the declared output')
            return asOutputType(((Map)value).get(outputName), outputType, workDir)
        }
        return asOutputType(requireJsonObject(value, 'record'), outputType, workDir)
    }

    /**
     * The final answer of an agent run, as returned by {@code stdout()} in an agent output.
     */
    static String answer(Object stdout) {
        return terminalAnswer(stdout)
    }

    /**
     * The {@code output} string of the {@code complete} frame a canonical task prints last on stdout.
     * Blank stdout is one failure; a last frame that is not an object, is not {@code complete}, or
     * carries no {@code output} is the other -- all three yield the same message, because from the
     * caller's side they are the same thing: no answer came back.
     */
    private static String terminalAnswer(Object stdout) {
        final String text = stdout?.toString()?.trim()
        if( !text )
            throw new ScriptRuntimeException('Canonical agent task completed without a result frame on stdout')
        final List<String> lines = text.readLines().findAll { it?.trim() }
        final Object frame = new JsonSlurper().parseText(lines.last())
        final Object answer = frame instanceof Map && ((Map)frame).get('type') == 'complete'
                ? ((Map)frame).get('output')
                : null
        if( answer == null )
            throw new ScriptRuntimeException('Canonical agent task returned an invalid terminal result frame')
        return answer.toString()
    }

    /** The decoded answer as a JSON object, or the {@code kind}-specific failure. */
    private static Map requireJsonObject(Object value, String kind) {
        if( !(value instanceof Map) )
            throw new ScriptRuntimeException("Canonical agent ${kind} output must be a JSON object")
        return (Map) value
    }

    /**
     * Bind an in-JVM runner's result into the task context.
     */
    void bind(Map ctx, Object result, AgentOutput output, Path workDir = null) {
        // no model-answered output: the agent's result is the files it wrote
        // and its final text, which is available through `stdout()`
        if( output == null )
            return
        final map = new JsonSlurper().parseText(stripFences(result as String)) as Map
        final value = mode == AgentOutputMode.SCALAR_CONTRACT ? map[output.name] : map
        ctx.put(output.name, asOutputType(value, output.type, workDir))
    }

    /**
     * Convert a decoded answer to the declared output type, resolving
     * relative paths against the agent work dir, since the model refers
     * to files by their name in the work dir.
     *
     * @param value
     * @param type
     * @param workDir
     */
    private static Object asOutputType(Object value, Type type, Path workDir) {
        return TypeHelper.asType(resolvePaths(value, type, workDir), type)
    }

    private static Object resolvePaths(Object value, Type type, Path workDir) {
        if( value == null )
            return value
        if( type == Path && value instanceof CharSequence ) {
            final str = value.toString()
            final path = workDir == null || str.contains('://') || Path.of(str).isAbsolute()
                ? FileHelper.asPath(str)
                : workDir.resolve(str)
            if( !Files.exists(path) )
                throw new ScriptRuntimeException("Agent output path '${str}' does not exist")
            return path
        }
        if( TypeHelper.isCollectionType(type) && value instanceof Collection ) {
            final elementType = type instanceof ParameterizedType ? ((ParameterizedType) type).getActualTypeArguments()[0] : Object
            return ((Collection) value).collect { el -> resolvePaths(el, elementType, workDir) }
        }
        if( TypeHelper.isRecordType(type) && value instanceof Map ) {
            final result = new LinkedHashMap((Map) value)
            for( final field : ((Class) type).getDeclaredFields() ) {
                if( !field.isSynthetic() && result.containsKey(field.getName()) )
                    result.put(field.getName(), resolvePaths(result.get(field.getName()), field.getGenericType(), workDir))
            }
            return result
        }
        return value
    }

    /**
     * Strip a leading ```json (or ```) fence and a trailing ``` fence from the
     * given text, returning the inner content. If no fences are present the
     * original text is returned unchanged.
     */
    private static String stripFences(String text) {
        if( text == null )
            return null
        final String s = text.trim()
        if( !s.startsWith('```') )
            return text
        // the opening fence line (``` or ```json) ends at the first newline; without one there is
        // no fenced block to strip
        final int open = s.indexOf('\n')
        if( open < 0 )
            return text
        // the LAST closing fence, so a fenced block that itself contains ``` keeps it; an
        // unterminated fence (no closing run at all) keeps everything after the opening line
        final int close = s.lastIndexOf('```')
        return (close > open ? s.substring(open + 1, close) : s.substring(open + 1)).trim()
    }
}
