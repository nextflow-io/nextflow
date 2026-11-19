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

import java.nio.file.Files

import nextflow.cli.project.CmdProjectList
import nextflow.cli.project.CmdProjectRemove
import nextflow.cli.project.CmdProjectView
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import spock.lang.Specification
import test.OutputCapture

class CmdProjectTest extends Specification {

    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    File root
    File savedRoot

    def setup() {
        savedRoot = AssetManager.root
        root = Files.createTempDirectory('test').toFile()
        AssetManager.root = root
    }

    def cleanup() {
        Plugins.stop()
        AssetManager.root = savedRoot
        root.deleteDir()
    }

    def 'should list the project subcommands'() {
        when:
        new CmdProject().usage()
        then:
        final out = capture.toString()
        out.contains('Usage: nextflow project <command> [options]')
        ['info', 'list', 'pull', 'remove', 'view'].every { name -> out.contains("  ${name.padRight(12)}") }
    }

    def 'should list the projects'() {
        when:
        new Launcher().parseMainArgs(args as String[]).command.run()
        then:
        capture.toString().contains('(none)')

        where:
        args << [['project', 'list'], ['list']]
    }

    def 'should run the project commands and their deprecated aliases'() {
        when:
        new Launcher().parseMainArgs(args as String[]).command.run()
        then:
        def e = thrown(AbortOperationException)
        e.message.startsWith(MESSAGE)

        where:
        args                                        | MESSAGE
        ['project', 'remove', 'acme/hello', '-f']   | 'No match found for: acme/hello'
        ['drop', 'acme/hello']                      | 'No match found for: acme/hello'
        ['project', 'info', 'acme/hello']           | 'Unknown project `acme/hello`'
        ['info', 'acme/hello']                      | 'Unknown project `acme/hello`'
        ['project', 'view', 'acme/hello']           | 'Unknown project `acme/hello'
        ['view', 'acme/hello']                      | 'Unknown project `acme/hello'
    }

    def 'should parse the deprecated top-level commands as aliases'() {
        when:
        def cmd = new Launcher().parseMainArgs('drop', 'acme/hello', '-r', 'v1', '-f').command
        then:
        cmd instanceof CmdDrop
        cmd instanceof CmdProjectRemove
        cmd.name == 'drop'
        (cmd as CmdDrop).args == ['acme/hello']
        (cmd as CmdDrop).revision == 'v1'
        (cmd as CmdDrop).force

        when:
        cmd = new Launcher().parseMainArgs('list').command
        then:
        cmd instanceof CmdList
        cmd instanceof CmdProjectList

        when:
        cmd = new Launcher().parseMainArgs('view', '-l', '-q', 'acme/hello').command
        then:
        cmd instanceof CmdView
        cmd instanceof CmdProjectView
        (cmd as CmdView).args == ['acme/hello']
        (cmd as CmdView).all
        (cmd as CmdView).quiet
    }

    def 'should hide the deprecated top-level commands'() {
        when:
        new Launcher().printCommands([new CmdDrop(), new CmdInfo(), new CmdList(), new CmdProject(), new CmdPull(), new CmdView()])
        then:
        capture.toString() == '''
                Commands:
                  info      Print system runtime information
                  project   Manage downloaded projects
                  pull      Download or update a project

                '''
                .stripIndent()
    }

    def 'should point to the new command in the help of a deprecated command'() {
        when:
        new Launcher().usage(NAME)
        then:
        capture.toString().contains("(deprecated -- use `nextflow project ${REPLACEMENT}`)")

        where:
        NAME    | REPLACEMENT
        'drop'  | 'remove'
        'list'  | 'list'
        'view'  | 'view'
    }
}
