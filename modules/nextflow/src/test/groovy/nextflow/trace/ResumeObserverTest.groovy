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

package nextflow.trace

import java.time.Clock
import java.time.Instant
import java.time.ZoneOffset

import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortSignalException
import spock.lang.Specification
import spock.lang.Unroll
import sun.misc.Signal

class ResumeObserverTest extends Specification {

    static final String TOKEN = 'dummy-token-replace-me'

    static final long PID = 74877

    static final Clock CLOCK = Clock.fixed(Instant.ofEpochSecond(1795543200), ZoneOffset.UTC)

    static final List<String> COMMAND = ['/opt/bin/nextflow', 'run', '-resume', '11111111-2222-3333-4444-555555555555', 'main.nf']


    private ResumeObserver observer(OutputStream out) {
        new ResumeObserver(TOKEN, out, PID, CLOCK, '/work/project')
    }

    private static AbortSignalException abort(String name) {
        new AbortSignalException(Signal.getConstructor(String).newInstance(name))
    }

    private Session sessionWith(List<String> command, Throwable error = null) {
        return Mock(Session) {
            getResumeCommand() >> command
            getError() >> error
        }
    }

    @Unroll
    def 'should be enabled=#expected with NXF_OSC7577=#opt, tty=#tty and token=#token' () {
        given:
        final env = [:]
        if( opt != null ) env.NXF_OSC7577 = opt
        if( token != null ) env.OSC7577_TOKEN = token
        SysEnv.push(env)

        expect:
        ResumeObserver.enabled(tty) == expected

        cleanup:
        SysEnv.pop()

        where:
        opt     | tty   | token   | expected
        null    | true  | 'abc'   | false
        'false' | true  | 'abc'   | false
        'true'  | false | 'abc'   | false
        'true'  | true  | null    | false
        'true'  | true  | ''      | false
        'true'  | true  | 'abc'   | true
        '1'     | true  | 'abc'   | true
    }

    def 'should report the resume command when the run is created' () {
        given:
        final out = new ByteArrayOutputStream()

        when:
        observer(out).onFlowCreate(sessionWith(COMMAND))

        then:
        out.toString('UTF-8') == Osc7577.report(TOKEN, PID, 1795543200, COMMAND, '/work/project')
        and:
        out.toString('UTF-8').contains(':cwd=' + Base64.encoder.encodeToString('/work/project'.bytes))
    }

    def 'should send nothing when there is no resume command' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(null))
        obs.onFlowComplete()

        then:
        out.size() == 0
    }

    def 'should send nothing when the resume command cannot be reported' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(['/opt/bin/nextflow', 'run', 'bad\narg']))
        obs.onFlowComplete()

        then:
        out.size() == 0
    }

    def 'should clear the record when the run completes' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(COMMAND))
        obs.onFlowComplete()

        then:
        out.toString('UTF-8') ==
            Osc7577.report(TOKEN, PID, 1795543200, COMMAND, '/work/project') +
            Osc7577.clear(TOKEN, PID, 1795543200)
    }

    def 'should clear the record when the run fails' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(COMMAND, new IllegalStateException('boom')))
        obs.onFlowComplete()

        then:
        out.toString('UTF-8').endsWith(Osc7577.clear(TOKEN, PID, 1795543200))
    }

    def 'should clear the record when the user interrupts with ctrl-c' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(COMMAND, abort('INT')))
        obs.onFlowComplete()

        then:
        out.toString('UTF-8').endsWith(Osc7577.clear(TOKEN, PID, 1795543200))
    }

    @Unroll
    def 'should keep the record when the run ends with SIG#signal' () {
        given:
        final out = new ByteArrayOutputStream()
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(COMMAND, abort(signal)))
        obs.onFlowComplete()

        then:
        out.toString('UTF-8') == Osc7577.report(TOKEN, PID, 1795543200, COMMAND, '/work/project')

        where:
        signal << ['HUP', 'TERM']
    }

    def 'should not fail the run when the terminal is gone' () {
        given:
        final out = new OutputStream() {
            void write(int b) { throw new IOException('Input/output error') }
        }
        final obs = observer(out)

        when:
        obs.onFlowCreate(sessionWith(COMMAND))
        obs.onFlowComplete()

        then:
        noExceptionThrown()
    }

}
