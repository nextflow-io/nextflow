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

import java.nio.file.Files
import java.nio.file.Path

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

    def 'should reject a scalar contract whose object lacks the declared output'() {
        when:
        plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame(output), 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent scalar output must be a JSON object containing the declared output'

        where:
        output << ['{"different":1}', '"a bare string"', '[1,2]', '42']
    }

    def 'should reject a record answer that is not a JSON object'() {
        when:
        plan(AgentOutputMode.RECORD).decode(frame('"a bare string"'), 'answer', String)

        then:
        final e = thrown(ScriptRuntimeException)
        e.message == 'Canonical agent record output must be a JSON object'
    }

    // --- paths ---------------------------------------------------------------------------------

    static class Report implements nextflow.script.types.Record {
        Path summary
        List<Path> files
    }

    static class Holder {
        public List<Path> files
    }

    def 'should resolve relative paths against the agent work dir'() {
        given:
        def workDir = Files.createTempDirectory('test')
        def other = Files.createTempFile('other', '.txt')
        Files.createFile(workDir.resolve('a.txt'))
        Files.createFile(workDir.resolve('b.txt'))

        when:
        final path = plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame('{"report":"a.txt"}'), 'report', Path, workDir)
        then:
        path == workDir.resolve('a.txt')

        when:
        final json = groovy.json.JsonOutput.toJson([summary: 'a.txt', files: ['b.txt', other.toString()]])
        final rec = plan(AgentOutputMode.RECORD).decode(frame(json), 'report', Report, workDir)
        then:
        rec.summary == workDir.resolve('a.txt')
        rec.files == [workDir.resolve('b.txt'), other]

        cleanup:
        workDir?.deleteDir()
        Files.deleteIfExists(other)
    }

    def 'should decode a list output with its element type'() {
        given:
        def workDir = Files.createTempDirectory('test')
        Files.createFile(workDir.resolve('a.txt'))
        final type = Holder.getField('files').getGenericType()

        when:
        final files = plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame('{"files":["a.txt"]}'), 'files', type, workDir)
        then:
        files == [workDir.resolve('a.txt')]

        cleanup:
        workDir?.deleteDir()
    }

    def 'should report a missing output path'() {
        given:
        def workDir = Files.createTempDirectory('test')

        when:
        plan(AgentOutputMode.SCALAR_CONTRACT).decode(frame('{"report":"missing.txt"}'), 'report', Path, workDir)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == "Agent output path 'missing.txt' does not exist"

        cleanup:
        workDir?.deleteDir()
    }
}
