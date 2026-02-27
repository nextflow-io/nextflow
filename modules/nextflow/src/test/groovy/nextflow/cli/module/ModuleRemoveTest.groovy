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

import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleStorage
import nextflow.pipeline.PipelineSpec
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for ModuleRemove command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleRemoveTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    // No setup needed - using root field directly

    def 'should remove module files and config entry'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = createTestModule(storage, reference)

        // Create spec file with module entry
        def specFile = new PipelineSpec(tempDir)
        specFile.addModuleEntry('@nf-core/fastqc', '1.0.0')

        and:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/fastqc']
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Removing module files')
        output.contains('Module files removed successfully')
        output.contains('Removing module entry from nextflow_spec.json')
        output.contains('Module entry removed from configuration')
        !Files.exists(moduleDir)

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@nf-core/fastqc') == null
    }

    def 'should keep config with -keep-config flag'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = createTestModule(storage, reference)

        // Create spec file
        def specFile = new PipelineSpec(tempDir)
        specFile.addModuleEntry('@nf-core/fastqc', '1.0.0')

        and:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/fastqc']
        cmd.keepConfig = true
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Removing module files')
        output.contains('Keeping module entry in nextflow_spec.json')
        !Files.exists(moduleDir)

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@nf-core/fastqc') == '1.0.0'
    }

    def 'should keep files with -keep-files flag'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = createTestModule(storage, reference)

        // Create spec file
        def specFile = new PipelineSpec(tempDir)
        specFile.addModuleEntry('@nf-core/fastqc', '1.0.0')

        and:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/fastqc']
        cmd.keepFiles = true
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Keeping module files')
        output.contains('Removing module entry from nextflow_spec.json')
        Files.exists(moduleDir)
        Files.exists(moduleDir.resolve('main.nf'))

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@nf-core/fastqc') == null
    }

    def 'should fail when both keep flags are set'() {
        given:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/fastqc']
        cmd.keepConfig = true
        cmd.keepFiles = true
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Cannot use both -keep-config and -keep-files')
    }

    def 'should handle removing non-existent module'() {
        given:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/nonexistent']
        cmd.root = tempDir

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('was not installed locally') || output.contains('was not found')
    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new ModuleRemove()
        cmd.args = []
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail with too many arguments'() {
        given:
        def cmd = new ModuleRemove()
        cmd.args = ['nf-core/fastqc', 'extra-arg']
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    private Path createTestModule(ModuleStorage storage, ModuleReference reference) {
        def moduleDir = storage.getModuleDir(reference)
        Files.createDirectories(moduleDir)

        // Create main.nf
        moduleDir.resolve('main.nf').text = '''
            process FASTQC {
                input:
                path reads

                output:
                path "*.html"

                script:
                """
                fastqc ${reads}
                """
            }
        '''.stripIndent()

        // Create meta.yml
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
            description: FastQC quality control
        '''.stripIndent()

        return moduleDir
    }
}
