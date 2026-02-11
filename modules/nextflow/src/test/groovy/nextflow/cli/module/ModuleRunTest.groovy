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

import io.seqera.npr.api.schema.v1.Module
import io.seqera.npr.api.schema.v1.ModuleRelease
import nextflow.cli.CliOptions
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleRegistryClient
import nextflow.module.ModuleStorage
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPOutputStream

/**
 * Tests for ModuleRun command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleRunTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    def 'should run module and create output file'() {
        given:
        // Create a simple module script that creates a file
        def moduleScript = '''
            process CREATE_FILE {
                output:
                path "test_output.txt"

                script:
                """
                echo "Module executed successfully" > test_output.txt
                """
            }
        '''.stripIndent()

        and:
        // Create module directory structure
        def storage = new ModuleStorage(tempDir)
        def moduleRef = new ModuleReference('nf-core', 'test-module')
        def moduleDir = storage.getModuleDir(moduleRef)
        Files.createDirectories(moduleDir)

        // Write main.nf
        moduleDir.resolve('main.nf').text = moduleScript

        // Write meta.yml
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/test-module
            version: 1.0.0
            description: Test module that creates a file
        '''.stripIndent()

        and:
        // Create mock module package
        def modulePackage = createModulePackage(moduleScript)

        // Mock registry client
        def mockClient = Stub(ModuleRegistryClient)
        def moduleRelease = new ModuleRelease()
        moduleRelease.version = '1.0.0'
        def module = new Module()
        module.name = '@nf-core/test-module'
        module.latest = moduleRelease
        mockClient.fetchModule(_) >> module  // Use wildcard to match any argument
        mockClient.downloadModule(_, _, _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }

        and:
        def cmd = new ModuleRun()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> new CliOptions()
            getCliString() >> "nextflow module run nf-core/test-module"
        }
        cmd.args = ['nf-core/test-module']
        cmd.root = tempDir
        cmd.client = mockClient

        when:
        cmd.run()
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }.join(" ")

        then:
        stdout.contains('Executing module...')
        stdout.contains('Process CREATE_FILE Outputs:')
        stdout.contains("test_output.txt")
        and:
        // Verify module was installed
        Files.exists(moduleDir)
        Files.exists(moduleDir.resolve('main.nf'))

    }

    def 'should run module with specific version'() {
        given:
        def moduleScript = '''
            process CREATE_FILE_V2 {
                output:
                path "test_output_v2.txt"

                script:
                """
                echo "Module version 2.0.0 executed successfully" > test_output_v2.txt
                """
            }
        '''

        and:
        def storage = new ModuleStorage(tempDir)
        def moduleRef = new ModuleReference('nf-core', 'test-module')
        def moduleDir = storage.getModuleDir(moduleRef)
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = moduleScript
        moduleDir.resolve('meta.yml').text = 'name: nf-core/test-module\nversion: 2.0.0'

        and:
        def modulePackage = createModulePackage(moduleScript)

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.downloadModule('@nf-core/test-module', '2.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }

        and:
        def cmd = new ModuleRun()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> new CliOptions()
            getCliString() >> "nextflow module run nf-core/test-module"
        }
        cmd.args = ['nf-core/test-module']
        cmd.version = '2.0.0'
        cmd.root = tempDir
        cmd.client = mockClient

        when:
        cmd.run()

        then:
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }.join(" ")
        stdout.contains('Executing module...')
        stdout.contains('Process CREATE_FILE_V2 Outputs:')
        stdout.contains("test_output_v2.txt")

    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new ModuleRun()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = []
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail with invalid module reference'() {
        given:
        def cmd = new ModuleRun()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['invalid-module']  // Missing scope
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    // Helper method to create a module package (tar.gz)
    private byte[] createModulePackage(String mainNfContent) {
        def baos = new ByteArrayOutputStream()

        new GZIPOutputStream(baos).withCloseable { gzos ->
            new TarArchiveOutputStream(gzos).withCloseable { tos ->
                // Add main.nf
                addTarEntry(tos, 'main.nf', mainNfContent.bytes)

                // Add meta.yml
                def metaContent = '''
                    name: test-module
                    version: 1.0.0
                    description: Test module
                '''.stripIndent()
                addTarEntry(tos, 'meta.yml', metaContent.bytes)
            }
        }

        return baos.toByteArray()
    }

    private void addTarEntry(TarArchiveOutputStream tos, String name, byte[] content) {
        def entry = new TarArchiveEntry(name)
        entry.setSize(content.length)
        tos.putArchiveEntry(entry)
        tos.write(content)
        tos.closeArchiveEntry()
    }
}
