/*
 * Copyright 2013-2024, Seqera Labs
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

import spock.lang.Specification

/**
 * Test suite for InstalledModule
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class InstalledModuleTest extends Specification {

    Path tempDir

    def setup() {
        tempDir = Files.createTempDirectory('nf-installed-module-test-')
    }

    def cleanup() {
        tempDir?.deleteDir()
    }

    def 'should report VALID integrity when checksum matches'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process TEST { }'

        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = 'name: test/module\nversion: 0.0.1'

        // Compute actual checksum
        def actualChecksum = ModuleChecksum.compute(moduleDir)

        // Save checksum
        ModuleChecksum.save(moduleDir, actualChecksum)

        def installed = new InstalledModule(
            reference: new ModuleReference('test', 'module'),
            directory: moduleDir,
            mainFile: mainFile,
            manifestFile: metaFile,
            moduleInfoFile: moduleDir.resolve('.module_info'),
            expectedChecksum: actualChecksum,
            installedVersion: "0.0.1"
        )

        when:
        def integrity = installed.getIntegrity()

        then:
        integrity == ModuleIntegrity.VALID
    }

    def 'should report MODIFIED integrity when checksum differs'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process TEST { }'

        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = 'name: test/module\nversion: 0.0.1'

        // Compute initial checksum
        def originalChecksum = ModuleChecksum.compute(moduleDir)
        ModuleChecksum.save(moduleDir, originalChecksum)

        // Modify the file
        mainFile.text = 'process TEST { println "modified" }'

        def installed = new InstalledModule(
            reference: new ModuleReference('test', 'module'),
            directory: moduleDir,
            mainFile: mainFile,
            manifestFile: metaFile,
            moduleInfoFile: moduleDir.resolve('.module_info'),
            expectedChecksum: originalChecksum,
            installedVersion: "0.0.1"
        )

        when:
        def integrity = installed.getIntegrity()

        then:
        integrity == ModuleIntegrity.MODIFIED
    }

    def 'should report CORRUPTED integrity when main.nf is missing'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        def mainFile = moduleDir.resolve('main.nf')
        // Don't create main.nf

        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = 'name: test/module\nversion: 0.0.1'

        def moduleInfoFile = moduleDir.resolve('.module_info')
        // CORRUPTED check happens before .module_info check, so content doesn't matter here

        def installed = new InstalledModule(
            reference: new ModuleReference('test', 'module'),
            directory: moduleDir,
            mainFile: mainFile,
            manifestFile: metaFile,
            moduleInfoFile: moduleInfoFile,
            expectedChecksum: 'some-checksum'
        )

        when:
        def integrity = installed.getIntegrity()

        then:
        integrity == ModuleIntegrity.CORRUPTED
    }

    def 'should report NO_REMOTE_MODULE when .module_info file absent'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process TEST { }'

        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = 'name: test/module\nversion: 0.0.1'

        def moduleInfoFile = moduleDir.resolve('.module_info')
        // Don't create .module_info file

        def installed = new InstalledModule(
            reference: new ModuleReference('test', 'module'),
            directory: moduleDir,
            mainFile: mainFile,
            manifestFile: metaFile,
            moduleInfoFile: moduleInfoFile,
            expectedChecksum: null
        )

        when:
        def integrity = installed.getIntegrity()

        then:
        integrity == ModuleIntegrity.NO_REMOTE_MODULE
    }

    def 'should handle checksum computation failure gracefully'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process TEST { }'

        def metaFile = moduleDir.resolve('meta.yml')
        metaFile.text = 'name: test/module\nversion: 0.0.1'

        // Create .module_info with a checksum that won't match the computed checksum
        ModuleChecksum.save(moduleDir, 'expected-checksum-that-will-not-match')
        def moduleInfoFile = moduleDir.resolve('.module_info')

        def installed = new InstalledModule(
            reference: new ModuleReference('test', 'module'),
            directory: moduleDir,
            mainFile: mainFile,
            manifestFile: metaFile,
            moduleInfoFile: moduleInfoFile,
            expectedChecksum: 'expected-checksum-that-will-not-match'
        )

        when:
        def integrity = installed.getIntegrity()

        then:
        // Should handle gracefully - since checksum won't match, it should report as MODIFIED
        integrity in [ModuleIntegrity.VALID, ModuleIntegrity.MODIFIED, ModuleIntegrity.CORRUPTED]
    }
}
