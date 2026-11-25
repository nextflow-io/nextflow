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
import spock.lang.Specification
import spock.lang.Unroll

class ProgramResumeTest extends Specification {

    static final String ID = '6f1c1b2e-8f3a-4c43-9d7e-0a8f4d2b9c11'

    def 'should resume a plain run on the same session'() {
        expect:
        ProgramResume.command('nextflow', ['run', 'main.nf'], ID) ==
                ['nextflow', 'run', '-resume', ID, 'main.nf']
    }

    def 'should keep top-level options before the run command'() {
        expect:
        ProgramResume.command('/usr/bin/nextflow', ['-c', 'site.config', '-q', 'run', 'main.nf', '-profile', 'test'], ID) ==
                ['/usr/bin/nextflow', '-c', 'site.config', '-q', 'run', '-resume', ID, 'main.nf', '-profile', 'test']
    }

    @Unroll
    def 'should replace the previous resume target: #name'() {
        expect:
        ProgramResume.command('nextflow', args, ID) == expected

        where:
        name                 | args                                                    | expected
        'resume last'        | ['run', 'main.nf', '-resume', 'last']                   | ['nextflow', 'run', '-resume', ID, 'main.nf']
        'resume other id'    | ['run', '-resume', 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeee', 'main.nf'] | ['nextflow', 'run', '-resume', ID, 'main.nf']
        'resume then params' | ['run', 'main.nf', '-resume', 'last', '--input', 'x']   | ['nextflow', 'run', '-resume', ID, 'main.nf', '--input', 'x']
    }

    def 'should drop the run name because a rerun cannot reuse it'() {
        expect:
        ProgramResume.command('nextflow', ['run', 'main.nf', '-name', 'happy_curie', '-with-trace', '-'], ID) ==
                ['nextflow', 'run', '-resume', ID, 'main.nf', '-with-trace', '-']
    }

    def 'should keep pipeline parameters that look like launcher options'() {
        expect:
        ProgramResume.command('nextflow', ['run', 'main.nf', '--name', 'sample', '--resume', 'x'], ID) ==
                ['nextflow', 'run', '-resume', ID, 'main.nf', '--name', 'sample', '--resume', 'x']
    }

    @Unroll
    def 'should not resume a command that is not a run: #args'() {
        expect:
        ProgramResume.command('nextflow', args, ID) == null

        where:
        args << [null, [], ['log'], ['pull', 'nf-core/rnaseq']]
    }

    @Unroll
    def 'should find the launcher from NXF_CLI #cli'() {
        expect:
        ProgramResume.launcher(cli, Paths.get('/home/me/project')) == expected

        where:
        cli                                    | expected
        null                                   | 'nextflow'
        ''                                     | 'nextflow'
        'nextflow run main.nf'                 | 'nextflow'
        '/usr/local/bin/nextflow run main.nf'  | '/usr/local/bin/nextflow'
        './nextflow run main.nf'               | '/home/me/project/nextflow'
        '../bin/nextflow run main.nf'          | '/home/me/bin/nextflow'
        '/opt/my tools/nextflow run main.nf'   | 'nextflow'
        '/opt/launch.sh run main.nf'           | 'nextflow'
    }

    @Unroll
    def 'should be off unless opted in: env=#env config=#config'() {
        given:
        SysEnv.push(env != null ? [NXF_OSC7577: env] : [:])
        def session = Stub(Session) {
            getConfig() >> (config != null ? [osc7577: config] : [:])
        }

        expect:
        ProgramResume.optedIn(session) == expected

        cleanup:
        SysEnv.pop()

        where:
        env     | config  | expected
        null    | null    | false
        null    | true    | true
        null    | false   | false
        null    | 'true'  | true
        'true'  | null    | true
        '1'     | null    | true
        'false' | null    | false
        'false' | true    | false
        'true'  | false   | true
        ''      | true    | true
        'yes'   | null    | false
    }

    def 'should withhold the token from a child only when opted in'() {
        given:
        SysEnv.push(env)
        def builder = new ProcessBuilder('true')
        builder.environment().put('OSC7577_TOKEN', 'dummy-token-replace-me')

        when:
        ProgramResume.withholdToken(builder, Stub(Session) { getConfig() >> [:] })

        then:
        builder.environment().containsKey('OSC7577_TOKEN') == !expected

        cleanup:
        SysEnv.pop()

        where:
        env                    | expected
        [NXF_OSC7577: 'true']  | true
        [:]                    | false
    }

}
