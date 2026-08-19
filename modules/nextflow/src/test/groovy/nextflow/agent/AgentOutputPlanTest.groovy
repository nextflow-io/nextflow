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

import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification

/**
 * {@link AgentOutputPlan#decode} and the fence stripping it depends on.
 *
 * <p>Both were previously untested: {@code decode} was only ever exercised end-to-end through a
 * running agent task, and {@code stripFences} not at all. Every branch below is one the MODEL can
 * put the driver into by answering in a slightly different shape, so the failure they guard is a
 * pipeline that dies decoding a correct answer.
 *
 * <p>{@code stripFences} is private and is deliberately tested THROUGH {@code decode} rather than
 * reflectively: what matters is that a fenced answer decodes, not how the fence is removed.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentOutputPlanTest extends Specification {

    /** The fence marker, spelled once so no formatter can reflow it inside a literal. */
    static private final String F = '```'

    static private String frame(String output) {
        return '{"type":"complete","output":' + groovy.json.JsonOutput.toJson(output) + '}'
    }

    static private AgentOutputPlan plan(AgentOutputMode mode) {
        return new AgentOutputPlan(mode, null)
    }

    // --- the terminal frame ------------------------------------------------------------------

    def 'should read the answer out of the terminal complete frame'() {
        expect:
        plan(AgentOutputMode.TEXT).decode(frame('hello'), 'answer', String) == 'hello'
    }

    def 'should read the LAST non-blank line, ignoring earlier frames and blank lines'() {
        given:
        final stdout = '{"type":"progress","step":1}\n\n' + frame('final') + '\n\n  \n'

        expect:
        plan(AgentOutputMode.TEXT).decode(stdout, 'answer', String) == 'final'
    }

    def 'should reject stdout that carries no result frame'() {
        when:
        plan(AgentOutputMode.TEXT).decode(stdout, 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent task completed without a result frame on stdout'

        where:
        stdout << [null, '', '   ', '\n\n']
    }

    def 'should reject a last frame that is not a complete result'() {
        when:
        plan(AgentOutputMode.TEXT).decode(stdout, 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent task returned an invalid terminal result frame'

        where:
        stdout << [
            '"just a string"',                          // not an object
            '[1,2,3]',                                  // not an object
            '{"type":"error","message":"boom"}',        // wrong type
            '{"type":"complete"}',                      // no output
            '{"type":"complete","output":null}',        // null output
            '{"type":"progress","output":"x"}',         // complete frame never arrived
        ]
    }

    // --- fence stripping, through decode -----------------------------------------------------

    def 'should decode a scalar answer whatever fence the model wrapped it in'() {
        expect:
        plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame(output), 'answer', Integer) == 42

        where:
        output << [
            '{"answer":42}',                            // no fence
            F + 'json\n{"answer":42}\n' + F,            // ```json fence
            F + '\n{"answer":42}\n' + F,                // bare fence
            F + 'json\n{"answer":42}',                  // UNTERMINATED fence
            '  ' + F + 'json\n{"answer":42}\n' + F + ' ',
        ]
    }

    def 'should leave a fence-less answer, and an unsplittable one-line fence, untouched'() {
        expect:
        // a ``` run with no newline after it is not a fenced block: there is nothing to strip, so
        // the text is handed to the JSON parser as-is (and fails there, not silently)
        plan(AgentOutputMode.TEXT).decode(frame(F + 'json {"answer":42} ' + F), 'answer', String) ==
                F + 'json {"answer":42} ' + F
    }

    // --- per-mode interpretation -------------------------------------------------------------

    def 'should unwrap the declared output for a scalar contract'() {
        expect:
        plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame('{"answer":"yes"}'), 'answer', String) == 'yes'
    }

    def 'should unwrap the declared output for a wrapped record'() {
        expect:
        plan(AgentOutputMode.WRAPPED).decode(frame('{"total":7,"other":1}'), 'total', Integer) == 7
    }

    def 'should reject a scalar contract whose object lacks the declared output'() {
        when:
        plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame(output), 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent scalar output must be a JSON object containing the declared output'

        where:
        output << ['{"different":1}', '"a bare string"', '[1,2]', '42']
    }

    def 'should reject a wrapped answer that is not a JSON object'() {
        when:
        plan(AgentOutputMode.WRAPPED).decode(frame('[1,2]'), 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent structured output must be a JSON object'
    }

    def 'should reject a record answer that is not a JSON object'() {
        when:
        plan(AgentOutputMode.RECORD).decode(frame('"a bare string"'), 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent record output must be a JSON object'
    }

    // --- mode predicates ---------------------------------------------------------------------

    def 'should report which modes are structured and which are wrapped'() {
        expect:
        plan(mode).isStructured() == structured
        plan(mode).isWrapped() == wrapped

        where:
        mode                              | structured | wrapped
        AgentOutputMode.TEXT              | false      | false
        AgentOutputMode.SCALAR_CONTRACT   | false      | false
        AgentOutputMode.RECORD            | true       | false
        AgentOutputMode.WRAPPED           | true       | true
    }
}
