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

import java.nio.file.Paths

import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortSignalException
import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
import spock.lang.Unroll
import sun.misc.Signal

class ProgramResumeObserverTest extends Specification {

    static final String TOKEN = 'dummy-token-replace-me'

    static final String ID = '6f1c1b2e-8f3a-4c43-9d7e-0a8f4d2b9c11'

    static final long PID = 4242

    static final long TIME = 1795543200

    static final Map<String,String> ENV = [OSC7577_TOKEN: TOKEN, NXF_CLI: '/usr/local/bin/nextflow run main.nf']

    static final List<String> ARGS = ['run', 'main.nf', '-resume', 'last', '-name', 'happy_curie']

    def setup() {
        SysEnv.push(ENV)
    }

    def cleanup() {
        SysEnv.pop()
    }

    private Session session(Map overrides = [:]) {
        final props = [cliArgs: ARGS, uniqueId: UUID.fromString(ID), preview: false, error: null] + overrides
        return Stub(Session) {
            getCliArgs() >> props.cliArgs
            getUniqueId() >> props.uniqueId
            isPreview() >> props.preview
            getError() >> props.error
        }
    }

    private ProgramResumeObserver observer(StringBuilder sink, boolean tty = true) {
        return new ProgramResumeObserver(sink, tty, PID, { TIME })
    }

    private static String cwd() {
        Paths.get('').toAbsolutePath().toString()
    }

    def 'should report how to resume the run when it starts'() {
        given:
        def sink = new StringBuilder()

        when:
        observer(sink).onFlowCreate(session())

        then:
        sink.toString() == Osc7577Emitter.report(PID, TIME,
                ['/usr/local/bin/nextflow', 'run', '-resume', ID, 'main.nf'],
                cwd(), TOKEN)
    }

    @Unroll
    def 'should report nothing when #why'() {
        given:
        SysEnv.push(ENV + env)
        def sink = new StringBuilder()

        when:
        def obs = observer(sink, tty)
        obs.onFlowCreate(session(sessionProps))
        obs.onFlowComplete()

        then:
        sink.toString() == ''

        cleanup:
        SysEnv.pop()

        where:
        why                              | tty   | env                          | sessionProps
        'stdout is not a terminal'       | false | [:]                          | [:]
        'the terminal gave no token'     | true  | [OSC7577_TOKEN: null]        | [:]
        'the token is empty'             | true  | [OSC7577_TOKEN: '']          | [:]
        'the token is malformed'         | true  | [OSC7577_TOKEN: 'too short'] | [:]
        'the run is a preview'           | true  | [:]                          | [preview: true]
        'there are no launcher arguments'| true  | [:]                          | [cliArgs: null]
        'the command is not a run'       | true  | [:]                          | [cliArgs: ['log']]
    }

    def 'should report nothing when the command does not fit the protocol'() {
        given:
        def sink = new StringBuilder()

        when:
        def obs = observer(sink)
        obs.onFlowCreate(session(cliArgs: ['run', 'main.nf', '--note', 'two\nlines']))
        obs.onFlowComplete()

        then:
        sink.toString() == ''
    }

    @Unroll
    def 'should clear the record when the run ends with #how'() {
        given:
        def sink = new StringBuilder()
        def obs = observer(sink)
        obs.onFlowCreate(session(error: error))
        sink.setLength(0)

        when:
        obs.onFlowComplete()

        then:
        sink.toString() == Osc7577Emitter.clear(PID, TIME, TOKEN)

        where:
        how                | error
        'success'          | null
        'a failure'        | new ScriptRuntimeException('boom')
        'ctrl-c'           | new AbortSignalException(new Signal('INT'))
    }

    @Unroll
    def 'should keep the record when the run ends with SIG#signal'() {
        given:
        def sink = new StringBuilder()
        def obs = observer(sink)
        obs.onFlowCreate(session(error: new AbortSignalException(new Signal(signal))))
        sink.setLength(0)

        when:
        obs.onFlowComplete()

        then:
        sink.toString() == ''

        where:
        signal << ['TERM', 'HUP']
    }

    def 'should expose the signal name without the SIG prefix'() {
        expect:
        new AbortSignalException(new Signal('HUP')).signalName == 'HUP'
    }

}
