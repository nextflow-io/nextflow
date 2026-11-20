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

package nextflow.cli.pipeline

import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.JCommander
import groovy.json.JsonSlurper
import nextflow.cli.CmdDrop
import nextflow.cli.CmdInfo
import nextflow.cli.CmdList
import nextflow.cli.CmdPipeline
import nextflow.cli.CmdPull
import nextflow.cli.CmdView
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

/**
 * Tests for the `-global` option of the pipeline commands and the deprecated project commands
 */
class CmdPipelineGlobalTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path base

    File savedRoot

    def setup() {
        savedRoot = AssetManager.root
        AssetManager.root = Files.createDirectories(base.resolve('assets')).toFile()
    }

    def cleanup() {
        Plugins.stop()
        AssetManager.root = savedRoot
    }

    private void seedGlobal(String... names) {
        names.each { Files.createDirectories(AssetManager.root.toPath().resolve(it)) }
    }

    // the printed lines, without the log
    private List<String> stdout() {
        return capture.toString().readLines().findAll { String it -> it && !it.contains(' DEBUG ') }
    }

    private static <T> T parse(String... args) {
        return new Launcher().parseMainArgs(args).command as T
    }

    def 'should install globally with pull, resolved against the primary registry'() {
        given:
        def install = new CmdPipelineInstall(launcher: new Launcher().parseMainArgs('pipeline'))
        new JCommander(install).parse('-g', 'nf-core/rnaseq', '-version', '3.14.0')

        when:
        def pull = install.pullCommand()
        then:
        install.global
        pull instanceof CmdPull
        pull.args == ['nf-core/rnaseq']
        pull.revision == '3.14.0'
        pull.hubProvider == 'registry.nextflow.io'
        pull.launcher == install.launcher
    }

    def 'should reject -force with a global install'() {
        when:
        new CmdPipelineInstall(args: ['nf-core/rnaseq'], global: true, force: true).run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Option `-force` is not supported with `-global`'
    }

    def 'should list the global pipelines'() {
        given:
        seedGlobal('nextflow-io/hello', 'acme/foo')

        when:
        new CmdPipelineList(global: true).run()
        then:
        stdout() == ['acme/foo', 'nextflow-io/hello']
    }

    def 'should list the global pipelines as json'() {
        given:
        seedGlobal('nextflow-io/hello', 'acme/foo')

        when:
        new CmdPipelineList(global: true, output: 'json').run()
        then:
        (new JsonSlurper().parseText(stdout().join('\n')) as Map).pipelines == [[name: 'acme/foo'], [name: 'nextflow-io/hello']]
    }

    def 'should report when no pipelines are installed globally'() {
        when:
        new CmdPipelineList(global: true).run()
        then:
        stdout() == ['No pipelines installed globally']
    }

    def 'should hint at the global pipelines when none are installed locally'() {
        given:
        seedGlobal(*GLOBAL)

        when:
        new CmdPipelineList(root: base.resolve('project')).run()
        then:
        stdout() == EXPECTED

        where:
        GLOBAL                              | EXPECTED
        []                                  | ['No pipelines installed']
        ['acme/foo']                        | ['No pipelines installed', '1 pipeline is installed globally -- use `nextflow pipeline list -g` to list them']
        ['acme/foo', 'nextflow-io/hello']   | ['No pipelines installed', '2 pipelines are installed globally -- use `nextflow pipeline list -g` to list them']
    }

    def 'should report a global pipeline that is not installed'() {
        when:
        new CmdPipelineRemove(args: ['acme/hello'], global: true).run()
        then:
        capture.toString().contains('Pipeline acme/hello is not installed globally')
    }

    def 'should view a global pipeline'() {
        when:
        new CmdPipelineView(args: ['acme/hello'], global: true).run()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Unknown project `acme/hello`'

        when:
        new CmdPipelineView(args: ['acme/hello'], global: true, version: '1.0').run()
        then:
        e = thrown(AbortOperationException)
        e.message == 'Option `-version` is not supported with `-global`'
    }

    def 'should hide the deprecated project commands'() {
        when:
        new Launcher().printCommands([new CmdDrop(), new CmdInfo(), new CmdList(), new CmdPipeline(), new CmdPull(), new CmdView()])
        then:
        capture.toString() == '''
                Commands:
                  info       Print system runtime information
                  pipeline   Manage registry pipelines
                  pull       Download or update a project

                '''
                .stripIndent()
    }

    def 'should point to the replacement in the help of a deprecated command'() {
        when:
        new Launcher().usage(NAME)
        then:
        capture.toString().contains(DESCRIPTION)

        where:
        NAME    | DESCRIPTION
        'drop'  | '(deprecated -- use `nextflow pipeline remove -g`)'
        'list'  | '(deprecated -- use `nextflow pipeline list -g`)'
        'view'  | '(deprecated)'
    }

    def 'should still run the deprecated project commands'() {
        given:
        seedGlobal('acme/foo')

        when:
        parse('list').run()
        then:
        capture.toString().contains('acme/foo')

        when:
        parse(*ARGS).run()
        then:
        def e = thrown(AbortOperationException)
        e.message.startsWith(MESSAGE)

        where:
        ARGS                        | MESSAGE
        ['drop', 'acme/hello']      | 'No match found for: acme/hello'
        ['info', 'acme/hello']      | 'Unknown project `acme/hello`'
        ['view', 'acme/hello']      | 'Unknown project `acme/hello'
    }
}
