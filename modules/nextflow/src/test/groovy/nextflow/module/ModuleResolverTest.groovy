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

import nextflow.config.ModulesConfig
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for ModuleResolver
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleResolverTest extends Specification {

    @TempDir
    Path tempDir

    def 'should throw exception when resolving non-installed module without auto-install'() {
        given:
        def resolver = new ModuleResolver(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')

        when:
        resolver.resolve(reference, null, false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('not installed')
        e.message.contains('nextflow module install')
    }

    def 'should throw exception when installed module is corrupted'() {
        given:
        def resolver = new ModuleResolver(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def storage = new ModuleStorage(tempDir)
        def moduleDir = storage.getModuleDir(reference)

        // Create corrupted module (directory exists but no main.nf)
        Files.createDirectories(moduleDir)
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
        '''

        when:
        resolver.resolve(reference, null, false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('corrupted')

        cleanup:
        FileHelper.deletePath(moduleDir)
    }

    def 'should warn about locally modified module'() {
        given:
        def resolver = new ModuleResolver(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def storage = new ModuleStorage(tempDir)
        def moduleDir = storage.getModuleDir(reference)

        // Create module with mismatched checksum
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
        '''
        moduleDir.resolve('.checksum').text = 'wrong-checksum'

        when:
        def result = resolver.resolve(reference, null, false)

        then:
        result != null
        result == moduleDir.resolve('main.nf')

        cleanup:
        FileHelper.deletePath(moduleDir)
    }

    def 'should throw exception when version mismatch without auto-install'() {
        given:
        def modulesConfig = new ModulesConfig(['@nf-core/fastqc': '2.0.0'])
        def resolver = new ModuleResolver(tempDir, modulesConfig, null)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def storage = new ModuleStorage(tempDir)
        def moduleDir = storage.getModuleDir(reference)

        // Create module with different version
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
        '''

        // Compute and save correct checksum
        def checksum = ModuleChecksum.compute(moduleDir)
        moduleDir.resolve('.checksum').text = checksum

        when:
        resolver.resolve(reference, null, false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('version mismatch')
        e.message.contains('installed=1.0.0')
        e.message.contains('required=2.0.0')

        cleanup:
        FileHelper.deletePath(moduleDir)
    }

    def 'should resolve installed module with matching version'() {
        given:
        def resolver = new ModuleResolver(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def storage = new ModuleStorage(tempDir)
        def moduleDir = storage.getModuleDir(reference)

        // Create valid module
        Files.createDirectories(moduleDir)
        def mainFile = moduleDir.resolve('main.nf')
        mainFile.text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
        '''

        // Compute and save correct checksum
        def checksum = ModuleChecksum.compute(moduleDir)
        moduleDir.resolve('.checksum').text = checksum

        when:
        def result = resolver.resolve(reference, '1.0.0', false)

        then:
        result == mainFile

        cleanup:
        FileHelper.deletePath(moduleDir)
    }

    def 'should throw exception when trying to update a module with local modifications without force'() {
        given:
        def resolver = new ModuleResolver(tempDir)
        def reference = new ModuleReference('nf-core', 'fastqc')
        def storage = new ModuleStorage(tempDir)
        def moduleDir = storage.getModuleDir(reference)

        // Create module with wrong checksum (simulating local modifications)
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = '''
            name: nf-core/fastqc
            version: 1.0.0
        '''
        moduleDir.resolve('.checksum').text = 'wrong-checksum'

        when:
        resolver.installModule(reference, '2.0.0', false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('local modifications')
        e.message.contains('--force')

        cleanup:
        FileHelper.deletePath(moduleDir)
    }
}
