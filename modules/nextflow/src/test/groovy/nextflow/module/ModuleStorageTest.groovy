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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPOutputStream
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream

import spock.lang.Specification

/**
 * Test suite for ModuleStorage
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleStorageTest extends Specification {

    Path tempDir

    def setup() {
        tempDir = Files.createTempDirectory('nf-module-storage-test-')
    }

    def cleanup() {
        tempDir?.deleteDir()
    }

    def 'should get module directory path'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')

        when:
        def moduleDir = storage.getModuleDir(reference)

        then:
        moduleDir == tempDir.resolve('modules/@nf-core/fastqc')
    }

    def 'should check if module is installed'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = storage.getModuleDir(reference)

        when:
        def installed = storage.isInstalled(reference)

        then:
        !installed

        when:
        Files.createDirectories(moduleDir)
        installed = storage.isInstalled(reference)

        then:
        installed
    }

    def 'should return null for non-existent module'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'nonexistent')

        when:
        def installed = storage.getInstalledModule(reference)

        then:
        installed == null
    }

    def 'should get installed module with metadata'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = storage.getModuleDir(reference)
        Files.createDirectories(moduleDir)

        // Create main.nf
        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process FASTQC { }'

        // Create meta.yml with version
        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = '''
            name: nf-core/fastqc
            version: 1.0.0
            description: FastQC quality control
            keywords:
              - quality-control
              - fastqc
        '''.stripIndent()

        // Create .checksum file
        def checksumFile = moduleDir.resolve('.checksum')
        checksumFile.text = 'abc123def456'

        when:
        def installed = storage.getInstalledModule(reference)

        then:
        installed != null
        installed.reference == reference
        installed.directory == moduleDir
        installed.mainFile == mainFile
        installed.manifestFile == moduleDir.resolve('meta.yml')
        installed.checksumFile == checksumFile
        installed.expectedChecksum == 'abc123def456'
        installed.installedVersion == '1.0.0'
    }

    def 'should list all installed modules'() {
        given:
        def storage = new ModuleStorage(tempDir)

        // Create multiple modules
        def modules = [
            new ModuleReference('nf-core', 'fastqc'),
            new ModuleReference('nf-core', 'multiqc'),
            new ModuleReference('myorg', 'custom')
        ]

        modules.each { ref ->
            def moduleDir = storage.getModuleDir(ref)
            Files.createDirectories(moduleDir)

            // Create main.nf
            moduleDir.resolve('main.nf').text = 'process TEST { }'

            // Create meta.yml with version
            moduleDir.resolve('meta.yml').text = """
                name: ${ref.nameWithoutPrefix}
                version: 1.0.0
            """.stripIndent()

            // Create .checksum
            moduleDir.resolve('.checksum').text = 'checksum'
        }

        when:
        def installed = storage.listInstalled()

        then:
        installed.size() == 3
        installed*.reference.fullName.sort() == ['@myorg/custom', '@nf-core/fastqc', '@nf-core/multiqc']
    }

    def 'should list nested modules recursively'() {
        given:
        def storage = new ModuleStorage(tempDir)

        // Create modules with nested paths
        def modules = [
            new ModuleReference('nf-core', 'fastqc'),
            new ModuleReference('nf-core', 'gfatools/gfa2fa'),
            new ModuleReference('nf-core', 'gfatools/gfa2gfa'),
            new ModuleReference('myorg', 'tools/subtools/module')
        ]

        modules.each { ref ->
            def moduleDir = storage.getModuleDir(ref)
            Files.createDirectories(moduleDir)

            // Create main.nf
            moduleDir.resolve('main.nf').text = 'process TEST { }'

            // Create meta.yml with version
            moduleDir.resolve('meta.yml').text = """
                name: ${ref.nameWithoutPrefix}
                version: 1.0.0
                description: Test module
                license: MIT
            """.stripIndent()

            // Create .checksum
            moduleDir.resolve('.checksum').text = 'checksum'
        }

        when:
        def installed = storage.listInstalled()

        then:
        installed.size() == 4
        installed*.reference.fullName.sort() == [
            '@myorg/tools/subtools/module',
            '@nf-core/fastqc',
            '@nf-core/gfatools/gfa2fa',
            '@nf-core/gfatools/gfa2gfa'
        ]
    }

    def 'should return empty list when no modules installed'() {
        given:
        def storage = new ModuleStorage(tempDir)

        when:
        def installed = storage.listInstalled()

        then:
        installed.isEmpty()
    }

    def 'should install module from gzip package'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def version = '1.0.0'

        // Create a gzip package file
        def packageFile = Files.createTempFile('module-', '.tgz')
        createTestPackage(packageFile)

        when:
        def installed = storage.installModule(reference, version, packageFile)

        then:
        installed != null
        installed.reference == reference
        installed.installedVersion == '1.0.0'
        Files.exists(installed.mainFile)
        Files.exists(installed.checksumFile)

        cleanup:
        packageFile?.delete()
    }

    def 'should replace existing module on install'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = storage.getModuleDir(reference)

        // Create existing installation
        Files.createDirectories(moduleDir)
        def oldFile = moduleDir.resolve('old-file.txt')
        oldFile.text = 'old content'

        // Create package
        def packageFile = Files.createTempFile('module-', '.tgz')
        createTestPackage(packageFile)

        when:
        def installed = storage.installModule(reference, '2.0.0', packageFile)

        then:
        installed != null
        !Files.exists(oldFile)  // Old file should be removed
        Files.exists(installed.mainFile)

        cleanup:
        packageFile?.delete()
    }

    def 'should remove installed module'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = storage.getModuleDir(reference)

        // Create module
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process TEST { }'

        expect:
        Files.exists(moduleDir)

        when:
        def removed = storage.removeModule(reference)

        then:
        removed
        !Files.exists(moduleDir)
    }

    def 'should return false when removing non-existent module'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'nonexistent')

        when:
        def removed = storage.removeModule(reference)

        then:
        !removed
    }

    def 'should compute and save checksum on install'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def packageFile = Files.createTempFile('module-', '.tgz')
        createTestPackage(packageFile)

        when:
        def installed = storage.installModule(reference, '1.0.0', packageFile)

        then:
        installed.expectedChecksum != null
        installed.expectedChecksum.length() > 0
        Files.exists(installed.checksumFile)

        cleanup:
        packageFile?.delete()
    }

    def 'should handle installation failure gracefully'() {
        given:
        def storage = new ModuleStorage(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def invalidPackage = Files.createTempFile('invalid-', '.tgz')
        invalidPackage.text = 'not a valid gzip file'

        when:
        storage.installModule(reference, '1.0.0', invalidPackage)

        then:
        thrown(Exception)

        and:
        !storage.isInstalled(reference)  // Should not leave partial installation

        cleanup:
        invalidPackage?.delete()
    }

    /**
     * Helper method to create a test package file
     */
    private void createTestPackage(Path packageFile) {
        // Create a temporary directory with module content
        def tempModuleDir = Files.createTempDirectory('temp-module-')

        // Create main.nf
        tempModuleDir.resolve('main.nf').text = '''
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
        tempModuleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
            description: FastQC quality control
            keywords:
              - quality-control
              - fastqc
        '''.stripIndent()

        // Create README
        tempModuleDir.resolve('README.md').text = '# FastQC Module'

        // Create tar.gz archive using Java libraries
        Files.newOutputStream(packageFile).withCloseable { fos ->
            new GZIPOutputStream(fos).withCloseable { gzos ->
                new TarArchiveOutputStream(gzos).withCloseable { tos ->
                    // Add all files from tempModuleDir
                    Files.walk(tempModuleDir).each { Path path ->
                        if (Files.isRegularFile(path)) {
                            // Get relative path
                            def relativePath = tempModuleDir.relativize(path).toString()

                            // Create tar entry
                            def entry = new TarArchiveEntry(path.toFile(), relativePath)
                            tos.putArchiveEntry(entry)

                            // Write file content
                            Files.copy(path, tos)

                            tos.closeArchiveEntry()
                        }
                    }
                }
            }
        }

        // Cleanup temp directory
        tempModuleDir.deleteDir()
    }
}
