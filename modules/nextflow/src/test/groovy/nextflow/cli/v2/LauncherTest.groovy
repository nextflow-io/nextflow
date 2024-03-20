/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cli.v2

import java.nio.file.Files

import picocli.CommandLine
import picocli.CommandLine.HelpCommand
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class LauncherTest extends Specification {

    def parseArgs(String... args) {
        def launcher = new Launcher()
        def cmd = new CommandLine(launcher)
        cmd.parseArgs(args)

        return launcher
    }

    def getCommand(Launcher launcher) {
        launcher.spec.commandLine().getParseResult().subcommand().commandSpec().commandLine().getCommand()
    }

    def getArgs(Launcher launcher) {
        launcher.spec.commandLine().getParseResult().originalArgs()
    }


    def 'should return `version` option' () {

        when:
        def launcher = parseArgs('-v')
        then:
        assert launcher.options.version

        when:
        launcher = parseArgs('--version')
        then:
        assert launcher.options.fullVersion

    }

    def 'should return `help` command' () {

        def launcher
        def command

        when:
        launcher = parseArgs('-h')
        then:
        assert launcher.help

        when:
        launcher = parseArgs('help')
        command = getCommand(launcher)
        then:
        command instanceof HelpCommand

        when:
        launcher = parseArgs('help','xxx')
        command = getCommand(launcher)
        then:
        command instanceof HelpCommand

    }

    def 'should return `info` command'() {

        when:
        def launcher = parseArgs('info')
        def command = getCommand(launcher)
        then:
        command instanceof InfoCmd
        command.pipeline == null

        when:
        launcher = parseArgs('info', 'xxx')
        command = getCommand(launcher)
        then:
        command instanceof InfoCmd
        command.pipeline == 'xxx'

    }

    def 'should return `pull` command'() {

        when:
        def launcher = parseArgs('pull', 'alpha')
        def command = getCommand(launcher)
        then:
        command instanceof PullCmd
        command.pipeline == 'alpha'

        when:
        launcher = parseArgs('pull', 'xxx', '--hub', 'bitbucket', '--user', 'xx:11')
        command = getCommand(launcher)
        then:
        command instanceof PullCmd
        command.pipeline == 'xxx'
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == '11'

    }

    def 'should return `clone` command'() {
        when:
        def launcher = parseArgs('clone', 'xxx', 'yyy')
        def command = getCommand(launcher)
        then:
        command instanceof CloneCmd
        command.pipeline == 'xxx'
        command.targetName == 'yyy'

        when:
        launcher = parseArgs('clone', 'xxx', '--hub', 'bitbucket', '--user', 'xx:yy')
        command = getCommand(launcher)
        then:
        command instanceof CloneCmd
        command.pipeline == 'xxx'
        command.targetName == null
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == 'yy'
    }

    def 'should return `run` command'() {
        when:
        def launcher = parseArgs('run', 'xxx', '--hub', 'bitbucket', '--user', 'xx:yy')
        def command = getCommand(launcher)
        then:
        command instanceof RunCmd
        command.pipeline == 'xxx'
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == 'yy'

        when:
        launcher = parseArgs('run', 'alpha', '--hub', 'github')
        command = getCommand(launcher)
        then:
        command instanceof RunCmd
        command.pipeline == 'alpha'
        command.hubProvider == 'github'

        when:
        launcher = parseArgs('run', 'script.nf', 'arg1', 'arg2', '--', '--alpha', '0', '--omega', '9')
        command = getCommand(launcher)
        then:
        command instanceof RunCmd
        command.pipeline == 'script.nf'
        command.args == ['arg1', 'arg2']
        command.params.'alpha' == '0'
        command.params.'omega' == '9'

    }

    def 'should make cli' () {
        given:
        def launcher = new Launcher()
        expect:
        launcher.makeCli('nf', 'run', 'foo.nf') == 'nf run foo.nf'
        launcher.makeCli('nf', 'run', 'foo.nf', '*.txt') == "nf run foo.nf '*.txt'"
        launcher.makeCli('/this/that/nf run foo.nf', 'run', 'foo.nf', 'a{1,2}.z') == "nf run foo.nf 'a{1,2}.z'"
        launcher.makeCli('/this/that/launch run bar.nf', 'run', 'bar.nf') == '/this/that/launch run bar.nf'
    }

}
