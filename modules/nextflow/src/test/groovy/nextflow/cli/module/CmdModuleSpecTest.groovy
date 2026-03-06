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

package nextflow.cli.module

import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for CmdModuleSpec CLI command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CmdModuleSpecTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    private CmdModuleSpec makeCmd(Path moduleDir = null, String scope = null) {
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        if( moduleDir )
            cmd.args = [tempDir.relativize(moduleDir).toString()]
        if( scope )
            cmd.moduleScope = scope
        return cmd
    }

    private Path createModule(String processCode, String dirName = 'my-module') {
        def moduleDir = Files.createDirectory(tempDir.resolve(dirName))
        def mainNf = moduleDir.resolve('main.nf')
        mainNf.text = processCode
        return moduleDir
    }

    /** Create a locally-installed module under tempDir/modules/@scope/name/ */
    private Path createInstalledModule(String scope, String name, String processCode) {
        def moduleDir = tempDir.resolve("modules/${scope}/${name}")
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = processCode
        return moduleDir
    }

    static final String SIMPLE_PROCESS = '''\
        process FASTQC {
            input:
            tuple val(meta), path(reads)

            output:
            path "*.html"

            script:
            "fastqc $reads"
        }
        '''.stripIndent()

    // =========================================================================
    // Argument resolution tests
    // =========================================================================

    def 'should throw when no args are provided'() {
        given:
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('required')
    }

    def 'should throw when path arg is used without -scope'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)   // no scope

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('-scope')
    }

    def 'should derive name from scope and process name when path arg is used with -scope'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir, 'nf-core')
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        capture.toString().contains('nf-core/fastqc')
    }

    def 'should resolve installed module reference and use scope/name as module name'() {
        given:
        createInstalledModule('nf-core', 'fastqc', SIMPLE_PROCESS)
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nf-core/fastqc']
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        capture.toString().contains('nf-core/fastqc')
    }

    def 'should ignore -scope when module reference arg is resolved'() {
        given:
        createInstalledModule('nf-core', 'fastqc', SIMPLE_PROCESS)
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nf-core/fastqc']
        cmd.moduleScope = 'other-scope'   // should be ignored
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('nf-core/fastqc')
        !output.contains('other-scope')
    }

    def 'should fall back to path resolution when module reference is not installed'() {
        given:
        // nf-core/fastqc is NOT in modules/ — arg is resolved as a relative path
        def moduleDir = tempDir.resolve('nf-core/fastqc')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = SIMPLE_PROCESS
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nf-core/fastqc']
        cmd.moduleScope = 'my-scope'
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        capture.toString().contains('my-scope/fastqc')
    }

    def 'should throw when module reference is not installed and no -scope is provided'() {
        given:
        // nf-core/fastqc is NOT installed; path nf-core/fastqc does not exist either
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nf-core/fastqc']

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    // =========================================================================
    // CLI tests
    // =========================================================================

    def 'should generate meta.yml for valid module path'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir, 'nf-core')

        when:
        cmd.run()

        then:
        def metaYml = moduleDir.resolve('meta.yml')
        Files.exists(metaYml)
        def content = metaYml.text
        content.contains('nf-core/fastqc')
        content.contains('input:')
        content.contains('output:')
        capture.toString().contains('Generated:')
    }

    def 'should generate meta.yml for installed module reference'() {
        given:
        def moduleDir = createInstalledModule('nf-core', 'fastqc', SIMPLE_PROCESS)
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nf-core/fastqc']

        when:
        cmd.run()

        then:
        def metaYml = moduleDir.resolve('meta.yml')
        Files.exists(metaYml)
        metaYml.text.contains('nf-core/fastqc')
        capture.toString().contains('Generated:')
    }

    def 'should throw AbortOperationException when main.nf is missing'() {
        given:
        def moduleDir = Files.createDirectory(tempDir.resolve('empty-module'))
        def cmd = makeCmd(moduleDir, 'test-scope')

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('main.nf')
    }

    def 'should throw AbortOperationException when directory does not exist'() {
        given:
        def cmd = new CmdModuleSpec()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        cmd.args = ['nonexistent-dir']

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should use first process and log warning when multiple processes are in main.nf'() {
        given:
        def multiProcessNf = '''\
            process FIRST {
                input:
                val x

                script:
                "echo $x"
            }

            process SECOND {
                input:
                path y

                script:
                "cat $y"
            }
            '''.stripIndent()
        def moduleDir = createModule(multiProcessNf)
        def cmd = makeCmd(moduleDir, 'test-scope')

        when:
        cmd.run()

        then:
        noExceptionThrown()
        def metaYml = moduleDir.resolve('meta.yml')
        Files.exists(metaYml)
        metaYml.text.contains('first')
    }

    def 'should throw AbortOperationException when no process in main.nf'() {
        given:
        def noProcessNf = '''\
            workflow {
                Channel.of(1) | view
            }
            '''.stripIndent()
        def moduleDir = createModule(noProcessNf)
        def cmd = makeCmd(moduleDir, 'test-scope')

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    // =========================================================================
    // Overwrite guard tests
    // =========================================================================

    def 'should throw AbortOperationException when meta.yml already exists and force is false'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def existingMetaYml = moduleDir.resolve('meta.yml')
        existingMetaYml.text = 'original content'
        def cmd = makeCmd(moduleDir, 'test-scope')
        cmd.force = false

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('-force')
        existingMetaYml.text == 'original content'
    }

    def 'should overwrite meta.yml when force flag is set'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def existingMetaYml = moduleDir.resolve('meta.yml')
        existingMetaYml.text = 'original content'
        def cmd = makeCmd(moduleDir, 'test-scope')
        cmd.force = true

        when:
        cmd.run()

        then:
        noExceptionThrown()
        existingMetaYml.text != 'original content'
        existingMetaYml.text.contains('fastqc')
    }

    // =========================================================================
    // Dry-run tests
    // =========================================================================

    def 'should print YAML to stdout and not write file when dry-run is set'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir, 'test-scope')
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        !Files.exists(moduleDir.resolve('meta.yml'))
        def output = capture.toString()
        output.contains('fastqc')
        output.contains('name:')
    }

    def 'should not modify existing meta.yml when dry-run is set'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def existingMetaYml = moduleDir.resolve('meta.yml')
        existingMetaYml.text = 'original content'
        def cmd = makeCmd(moduleDir, 'test-scope')
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        existingMetaYml.text == 'original content'
    }

    def 'should output correct YAML header when dry-run is used'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir, 'test-scope')
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('# This file was auto-generated')
        output.contains('name:')
    }

    def 'should use provided CLI optional fields in generated YAML'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir, 'nf-core')
        cmd.moduleVersion = '1.0.0'
        cmd.description = 'Run FastQC on reads'
        cmd.license = 'MIT'
        cmd.authors = ['Alice <alice@example.com>', 'Bob <bob@example.com>']
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('nf-core/fastqc')   // name derived from scope + process name
        output.contains('1.0.0')
        output.contains('Run FastQC on reads')
        output.contains('MIT')
        output.contains('Alice <alice@example.com>')
        output.contains('Bob <bob@example.com>')
    }
}
