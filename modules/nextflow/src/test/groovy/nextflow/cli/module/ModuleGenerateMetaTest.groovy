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
 * Tests for ModuleGenerateMeta CLI command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleGenerateMetaTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    private ModuleGenerateMeta makeCmd(Path moduleDir = null) {
        def cmd = new ModuleGenerateMeta()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir
        if( moduleDir ) {
            cmd.args = [tempDir.relativize(moduleDir).toString()]
        }
        return cmd
    }

    private Path createModule(String processCode, String dirName = 'my-module') {
        def moduleDir = Files.createDirectory(tempDir.resolve(dirName))
        def mainNf = moduleDir.resolve('main.nf')
        mainNf.text = processCode
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
    // T010: US1 CLI tests
    // =========================================================================

    def 'should generate meta.yml for valid module with args'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)

        when:
        cmd.run()

        then:
        def metaYml = moduleDir.resolve('meta.yml')
        Files.exists(metaYml)
        def content = metaYml.text
        content.contains('fastqc')
        content.contains('input:')
        content.contains('output:')
        capture.toString().contains('Generated:')
    }

    def 'should generate meta.yml in current dir when no args'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = SIMPLE_PROCESS
        def cmd = new ModuleGenerateMeta()
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        def metaYml = tempDir.resolve('meta.yml')
        Files.exists(metaYml)
    }

    def 'should throw AbortOperationException when main.nf is missing'() {
        given:
        def moduleDir = Files.createDirectory(tempDir.resolve('empty-module'))
        def cmd = makeCmd(moduleDir)

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('main.nf')
    }

    def 'should throw AbortOperationException when directory does not exist'() {
        given:
        def cmd = new ModuleGenerateMeta()
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
        def cmd = makeCmd(moduleDir)

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
        def cmd = makeCmd(moduleDir)

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    // =========================================================================
    // T012: US2 overwrite guard tests
    // =========================================================================

    def 'should throw AbortOperationException when meta.yml already exists and force is false'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def existingMetaYml = moduleDir.resolve('meta.yml')
        existingMetaYml.text = 'original content'
        def cmd = makeCmd(moduleDir)
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
        def cmd = makeCmd(moduleDir)
        cmd.force = true

        when:
        cmd.run()

        then:
        noExceptionThrown()
        existingMetaYml.text != 'original content'
        existingMetaYml.text.contains('fastqc')
    }

    // =========================================================================
    // T014: US3 dry-run tests
    // =========================================================================

    def 'should print YAML to stdout and not write file when dry-run is set'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)
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
        def cmd = makeCmd(moduleDir)
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        existingMetaYml.text == 'original content'
    }

    def 'should output correct YAML for module when dry-run is used'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('# This file was auto-generated')
        output.contains('name:')
    }

    def 'should use provided CLI args for required fields'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)
        cmd.moduleName = 'nf-core/fastqc'
        cmd.moduleVersion = '1.0.0'
        cmd.description = 'Run FastQC on reads'
        cmd.license = 'MIT'
        cmd.authors = ['Alice <alice@example.com>', 'Bob <bob@example.com>']
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('nf-core/fastqc')
        output.contains('1.0.0')
        output.contains('Run FastQC on reads')
        output.contains('MIT')
        output.contains('Alice <alice@example.com>')
        output.contains('Bob <bob@example.com>')
        // Required fields have no placeholders;
        !output.contains('TODO: Add version')
        !output.contains('TODO: Add module description')
        !output.contains('TODO: Add license')
        !output.contains('TODO: Add author')
    }

    def 'should fall back to TODO placeholders when CLI args are not provided'() {
        given:
        def moduleDir = createModule(SIMPLE_PROCESS)
        def cmd = makeCmd(moduleDir)
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('TODO: Add version')
        output.contains('TODO: Add module description')
        output.contains('TODO: Add license')
        output.contains('TODO: Add author')
    }
}
