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

import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Path

/**
 * Test CmdAuth functionality
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class CmdAuthTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    def 'should have correct name'() {
        given:
        def cmd = new CmdAuth()

        expect:
        cmd.getName() == 'auth'
    }


    def 'should show usage when no args provided'() {
        given:
        def cmd = new CmdAuth()

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('Manage Seqera Platform authentication')
        output.contains('Usage: nextflow auth <sub-command> [options]')
        output.contains('Commands:')
        output.contains('login')
        output.contains('logout')
        output.contains('config')
        output.contains('status')
    }

    def 'should show specific command usage'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['login']

        when:
        cmd.usage()

        then:
        def output = capture.toString()
        output.contains('Authenticate with Seqera Platform')
        output.contains('Usage: nextflow auth login')
        output.contains('-u, -url <endpoint>')
    }

    def 'should throw error for unknown command'() {
        given:
        def cmd = Spy(CmdAuth)
        cmd.args = ['unknown']
        def operation = Mock(CmdAuth.AuthCommand)

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        def ex = thrown(AbortOperationException)
        ex.message.contains('Unknown auth sub-command: unknown')
    }

    def 'should suggest closest command for typos'() {
        given:
        def cmd = new CmdAuth()

        when:
        cmd.getCmd(['loginn'])

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Unknown auth sub-command: loginn')
        ex.message.contains('Did you mean one of these?')
        ex.message.contains('login')
    }

    def 'login command should validate too many arguments'() {
        given:
        def cmd = Spy(CmdAuth)
        def operation = Mock(CmdAuth.AuthCommand)
        cmd.args = ['login', 'extra']

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for login command')
    }

    def 'logout command should validate too many arguments'() {
        given:
        def cmd = Spy(CmdAuth)
        def operation = Mock(CmdAuth.AuthCommand)
        cmd.args = ['logout', 'extra']

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for logout command')
    }

    def 'config command should validate too many arguments'() {
        given:
        def cmd = Spy(CmdAuth)
        def operation = Mock(CmdAuth.AuthCommand)
        cmd.args = ['config', 'extra']

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for config command')
    }

    def 'status command should validate too many arguments'() {
        given:
        def cmd = Spy(CmdAuth)
        def operation = Mock(CmdAuth.AuthCommand)
        cmd.args = ['status', 'extra']

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for status command')
    }

    def 'login command should use provided API URL'() {
        given:
        def cmd = Spy(CmdAuth)
        def operation = Mock(CmdAuth.AuthCommand)
        cmd.args = ['login']
        cmd.apiUrl = 'https://api.example.com'

        when:
        def loginCmd = cmd.getCmd(cmd.args)

        then:
        loginCmd instanceof CmdAuth.LoginCmd

        when:
        cmd.run()

        then:
        1 * cmd.loadOperation() >> operation
        loginCmd.apiUrl == 'https://api.example.com'
    }

    def 'should have all required subcommands'() {
        given:
        def cmd = new CmdAuth()

        expect:
        cmd.commands.size() == 4
        cmd.commands.find { it.name == 'login' } != null
        cmd.commands.find { it.name == 'logout' } != null
        cmd.commands.find { it.name == 'config' } != null
        cmd.commands.find { it.name == 'status' } != null
    }
}
