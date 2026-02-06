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
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleRegistryClient
import nextflow.pipeline.PipelineSpec
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
 * Tests for ModuleInstall command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleInstallTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    def 'should install module with latest version'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.root = tempDir

        and:
        // Create mock module package
        def modulePackage = createModulePackage('nf-core', 'fastqc', '1.0.0')

        // Mock registry client
        def mockClient = Mock(ModuleRegistryClient)
        mockClient.fetchModule('@nf-core/fastqc') >> new Module(
            name: '@nf-core/fastqc',
            latest: new ModuleRelease(version: '1.0.0')
        )
        mockClient.downloadModule('@nf-core/fastqc', '1.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Installing')
        output.contains('nf-core/fastqc')
        output.contains('1.0.0')

        and:
        def moduleDir = tempDir.resolve('modules/@nf-core/fastqc')
        Files.exists(moduleDir)
        Files.exists(moduleDir.resolve('main.nf'))
        Files.exists(moduleDir.resolve('meta.yml'))

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@nf-core/fastqc') == '1.0.0'
    }

    def 'should install module with specific version'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.version = '2.0.0'
        cmd.root = tempDir

        and:
        def modulePackage = createModulePackage('nf-core', 'fastqc', '2.0.0')

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.downloadModule(_, _, _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Installing')
        output.contains('nf-core/fastqc')
        output.contains('2.0.0')

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@nf-core/fastqc') == '2.0.0'

    }

    def 'should update existing module with force flag'() {
        given:
        // Pre-install version 1.0.0
        def moduleDir = tempDir.resolve('modules/@nf-core/fastqc')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process OLD { }'
        moduleDir.resolve('meta.yml').text = """
                                                    name: nf-core/fastqc
                                                    version: '1.0.0'
                                                    description: Test module
                                                    """.stripIndent()

        def spec = new PipelineSpec(tempDir)
        spec.addModuleEntry('@nf-core/fastqc', '1.0.0')

        and:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.version = '2.0.0'
        cmd.force = true
        cmd.root = tempDir

        and:
        def modulePackage = createModulePackage('nf-core', 'fastqc', '2.0.0')

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.downloadModule('@nf-core/fastqc', '2.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Installing')
        output.contains('2.0.0')

        and:
        def updatedSpec = new PipelineSpec(tempDir)
        updatedSpec.getModules().get('@nf-core/fastqc') == '2.0.0'

        and:
        moduleDir.resolve('main.nf').text.contains('FASTQC')  // New content
    }

    def 'should fail when module already installed without force'() {
        given:
        // Pre-install the module
        def moduleDir = tempDir.resolve('modules/@nf-core/fastqc')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process FASTQC { }'
        moduleDir.resolve('meta.yml').text = """
                                                    name: nf-core/fastqc
                                                    version: '1.0.0'
                                                    description: Test module
                                                    """.stripIndent()
        moduleDir.resolve('.checksum').text = 'wrong-checksum'
        def spec = new PipelineSpec(tempDir)
        spec.addModuleEntry('@nf-core/fastqc', '1.0.0')

        and:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.version = '2.0.0'
        cmd.root = tempDir

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.fetchModule('@nf-core/fastqc') >> new Module(
            name: '@nf-core/fastqc',
            latest: new ModuleRelease(version: '2.0.0')
        )
        cmd.client = mockClient

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('already installed') || e.message.contains('-force')
    }

    def 'should handle module with scope in name'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['myorg/custom-module']
        cmd.root = tempDir

        and:
        def modulePackage = createModulePackage('myorg', 'custom-module', '1.0.0')

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.fetchModule('@myorg/custom-module') >> new Module(
            name: '@myorg/custom-module',
            latest: new ModuleRelease(version: '1.0.0')
        )
        mockClient.downloadModule('@myorg/custom-module', '1.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()

        then:
        def moduleDir = tempDir.resolve('modules/@myorg/custom-module')
        Files.exists(moduleDir)
        Files.exists(moduleDir.resolve('main.nf'))

        and:
        def spec = new PipelineSpec(tempDir)
        spec.getModules().get('@myorg/custom-module') == '1.0.0'
    }

    def 'should create modules directory if it does not exist'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.root = tempDir

        and:
        def modulePackage = createModulePackage('nf-core', 'fastqc', '1.0.0')

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.fetchModule('@nf-core/fastqc') >> new Module(
            name: 'nf-core/fastqc',
            latest: new ModuleRelease(version: '1.0.0')
        )
        mockClient.downloadModule('@nf-core/fastqc', '1.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()

        then:
        Files.exists(tempDir.resolve('modules'))
        Files.exists(tempDir.resolve('modules/@nf-core'))
        Files.exists(tempDir.resolve('modules/@nf-core/fastqc'))
    }

    def 'should create checksum file after installation'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc']
        cmd.root = tempDir

        and:
        def modulePackage = createModulePackage('nf-core', 'fastqc', '1.0.0')

        def mockClient = Mock(ModuleRegistryClient)
        mockClient.fetchModule('@nf-core/fastqc') >> new Module(
            name: 'nf-core/fastqc',
            latest: new ModuleRelease(version: '1.0.0')
        )
        mockClient.downloadModule('@nf-core/fastqc', '1.0.0', _) >> { String name, String version, Path dest ->
            Files.write(dest, modulePackage)
            return dest
        }
        cmd.client = mockClient

        when:
        cmd.run()

        then:
        def moduleDir = tempDir.resolve('modules/@nf-core/fastqc')
        Files.exists(moduleDir.resolve('.checksum'))

        and:
        def checksum = moduleDir.resolve('.checksum').text
        checksum != null
        !checksum.isEmpty()
    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new ModuleInstall()
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

    def 'should fail with too many arguments'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['nf-core/fastqc', 'extra-arg']
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail with invalid module reference'() {
        given:
        def cmd = new ModuleInstall()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['invalid-module-name']  // Missing scope
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    // Helper method to create a module package (tar.gz)
    private byte[] createModulePackage(String scope, String name, String version) {
        def baos = new ByteArrayOutputStream()

        // Create tar.gz with module files
        def mainNfContent = """
            process ${name.toUpperCase().replaceAll('-', '_')} {
                input:
                path reads

                output:
                path "*.html"

                script:
                \"\"\"
                echo "Running ${name}"
                \"\"\"
            }
        """.stripIndent()
        new GZIPOutputStream(baos).withCloseable { gzos ->
            new TarArchiveOutputStream(gzos).withCloseable { tos ->
                // Add main.nf
                addTarEntry(tos, 'main.nf', mainNfContent.bytes)

                // Add meta.yml
                def metaContent = """
                    name: ${scope}/${name}
                    version: ${version}
                    description: Test module
                """.stripIndent()
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
