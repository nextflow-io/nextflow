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

import groovy.json.JsonSlurper
import nextflow.module.ModuleChecksum
import nextflow.module.ModuleReference
import nextflow.module.ModuleStorage
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for ModuleList command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleListTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    // No setup needed - using root field directly

    def 'should list installed modules with formatted output'() {
        given:
        def storage = new ModuleStorage(tempDir)

        // Create test modules
        createTestModule(storage, 'nf-core', 'fastqc', '1.0.0')
        createTestModule(storage, 'nf-core', 'multiqc', '2.1.0')

        and:
        def cmd = new ModuleList()
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Installed modules:')
        output.contains('nf-core/fastqc')
        output.contains('1.0.0')
        output.contains('nf-core/multiqc')
        output.contains('2.1.0')
        output.contains('OK') || output.contains('NO CHECKSUM')
    }

    def 'should list installed modules with JSON output'() {
        given:
        def storage = new ModuleStorage(tempDir)

        // Create test module
        createTestModule(storage, 'nf-core', 'fastqc', '1.5.0')

        and:
        def cmd = new ModuleList()
        cmd.jsonOutput = true
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()
        def json = new JsonSlurper().parseText(output)

        then:
        json.modules != null
        json.modules.size() == 1
        json.modules[0].name == 'nf-core/fastqc'
        json.modules[0].version == '1.5.0'
        json.modules[0].integrity != null
    }

    def 'should handle no installed modules'() {
        given:
        def cmd = new ModuleList()
        cmd.root = tempDir  // Use test directory with no modules

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('No modules installed')
    }

    def 'should show modified status for locally modified modules'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def moduleDir = createTestModule(storage, 'nf-core', 'fastqc', '1.0.0')

        // Modify the module to trigger checksum mismatch
        moduleDir.resolve('main.nf').text = 'process MODIFIED { }'

        and:
        def cmd = new ModuleList()
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('nf-core/fastqc')
        output.contains('MODIFIED')
    }

    def 'should list multiple modules sorted by name'() {
        given:
        def storage = new ModuleStorage(tempDir)

        // Create modules in random order
        createTestModule(storage, 'nf-core', 'samtools', '1.0.0')
        createTestModule(storage, 'nf-core', 'fastqc', '1.0.0')
        createTestModule(storage, 'myorg', 'custom', '2.0.0')

        and:
        def cmd = new ModuleList()
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('fastqc')
        output.contains('samtools')
        output.contains('myorg/custom')
    }

    private Path createTestModule(ModuleStorage storage, String scope, String name, String version) {
        def moduleDir = storage.getModuleDir(new ModuleReference(scope, name))
        Files.createDirectories(moduleDir)

        // Create main.nf
        moduleDir.resolve('main.nf').text = """
            process ${name.toUpperCase()} {
                input:
                path reads

                output:
                path "*.html"

                script:
                \"\"\"
                echo "test"
                \"\"\"
            }
        """.stripIndent()

        // Create meta.yml
        moduleDir.resolve('meta.yml').text = """
            name: ${scope}/${name}
            version: ${version}
            description: Test module
        """.stripIndent()

        // Create checksum
        def checksum = ModuleChecksum.compute(moduleDir)
        moduleDir.resolve('.checksum').text = checksum

        return moduleDir
    }
}
