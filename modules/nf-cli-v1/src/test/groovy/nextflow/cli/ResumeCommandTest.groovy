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

package nextflow.cli

import java.nio.file.Paths

import spock.lang.Specification
import spock.lang.Unroll

class ResumeCommandTest extends Specification {

    static final UUID ID = UUID.fromString('11111111-2222-3333-4444-555555555555')

    @Unroll
    def 'should build the resume command for #args' () {
        expect:
        ResumeCommand.of('/usr/bin/nextflow', args, ID) == expected

        where:
        args                                                        | expected
        ['run', 'main.nf']                                          | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'main.nf']
        ['run', 'nf-core/rnaseq', '-r', '3.14', '-profile', 'test'] | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'nf-core/rnaseq', '-r', '3.14', '-profile', 'test']
        ['run', 'main.nf', '-resume', 'last']                       | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'main.nf']
        ['run', 'main.nf', '-resume', 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee'] | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'main.nf']
        ['run', '-name', 'happy_darwin', 'main.nf']                 | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'main.nf']
        ['-q', '-c', 'my.config', 'run', 'main.nf']                 | ['/usr/bin/nextflow', '-q', '-c', 'my.config', 'run', '-resume', ID.toString(), 'main.nf']
        ['run', 'main.nf', '--input', 'my sample.csv', '--opt', '$(reboot)'] | ['/usr/bin/nextflow', 'run', '-resume', ID.toString(), 'main.nf', '--input', 'my sample.csv', '--opt', '$(reboot)']
    }

    def 'should build the resume command from the arguments the launcher normalizes' () {
        given:
        final launcher = new Launcher().parseMainArgs('run', 'main.nf', '-resume', '-name', 'my_run', '--foo', 'bar')
        expect:
        launcher.normalizedArgs == ['run', 'main.nf', '-resume', 'last', '-name', 'my_run', '--foo=bar']
        and:
        ResumeCommand.of('nextflow', launcher.normalizedArgs, ID) == ['nextflow', 'run', '-resume', ID.toString(), 'main.nf', '--foo=bar']
    }

    @Unroll
    def 'should not build a resume command for #args' () {
        expect:
        ResumeCommand.of('/usr/bin/nextflow', args, ID) == null

        where:
        args                                  | _
        []                                    | _
        ['info']                              | _
        ['module', 'list']                    | _
        ['run', '-']                          | _
        ['run', '-', '-with-trace']           | _
    }

    @Unroll
    def 'should find the launcher path for #cli' () {
        expect:
        ResumeCommand.launcherPath(cli, ['run', 'main.nf'] as String[], Paths.get('/work/project')) == expected

        where:
        cli                                | expected
        '/opt/bin/nextflow run main.nf'    | '/opt/bin/nextflow'
        'nextflow run main.nf'             | 'nextflow'
        './nextflow run main.nf'           | '/work/project/nextflow'
        '../tools/nextflow run main.nf'    | '/work/tools/nextflow'
        '/my tools/nextflow run main.nf'   | '/my tools/nextflow'
        null                               | 'nextflow'
        ''                                 | 'nextflow'
        'nextflow run other.nf'            | 'nextflow'
    }

    def 'should find the launcher path without arguments' () {
        expect:
        ResumeCommand.launcherPath('/opt/bin/nextflow ', [] as String[], Paths.get('/work')) == '/opt/bin/nextflow'
    }

}
