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

import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import nextflow.extension.Bolts
import org.slf4j.LoggerFactory
import nextflow.module.spi.RemoteModuleResolver
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
        ModuleChecksum.save(moduleDir, 'wrong-checksum')

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
        def resolver = new ModuleResolver(tempDir, new RegistryConfig())
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
        ModuleChecksum.save(moduleDir, checksum)

        when:
        resolver.resolve(reference, '2.0.0', false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('version mismatch')
        e.message.contains('installed=1.0.0')
        e.message.contains('required=2.0.0')

        cleanup:
        if( moduleDir ) FileHelper.deletePath(moduleDir)
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
        ModuleChecksum.save(moduleDir, checksum)

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
        ModuleChecksum.save(moduleDir, 'wrong-checksum')

        when:
        resolver.installWithDependencies(reference, '2.0.0', false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('local modifications')
        e.message.contains('-force')

        cleanup:
        FileHelper.deletePath(moduleDir)
    }

    def 'should read the pinned dependency version from the module spec'() {
        given:
        def moduleDir = tempDir.resolve('modules/nf-core/my-workflow')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('meta.yml').text = """\
            name: my-workflow
            requires:
              nextflow: '>=26.04.0'
              modules:
                - nf-core/fastqc@1.2.3
                - nf-core/multiqc
            """.stripIndent()

        expect:
        ModuleResolver.pinnedVersion(moduleDir, new ModuleReference('nf-core', 'fastqc')) == '1.2.3'
        and: 'a dependency declared without a version is not pinned'
        ModuleResolver.pinnedVersion(moduleDir, new ModuleReference('nf-core', 'multiqc')) == null
        and: 'a module that is not a declared dependency is not pinned'
        ModuleResolver.pinnedVersion(moduleDir, new ModuleReference('nf-core', 'bwa')) == null
        and: 'a directory without a module spec has no pins'
        ModuleResolver.pinnedVersion(tempDir, new ModuleReference('nf-core', 'fastqc')) == null
    }


    def 'does not warn for a module without .module-info inside a pipeline root'() {
        given:
        // the warn1 cache is JVM-global and would swallow the message
        def cache = Bolts.getDeclaredField('LOGGER_CACHE')
        cache.setAccessible(true)
        (cache.get(null) as Map).clear()
        def logger = (Logger) LoggerFactory.getLogger(ModuleResolver)
        def level = logger.level
        logger.setLevel(Level.WARN)
        def appender = new ListAppender<ILoggingEvent>()
        appender.start()
        logger.addAppender(appender)

        and:
        def reference = new ModuleReference('nf-core', 'fastqc')
        def moduleDir = new ModuleStorage(tempDir).getModuleDir(reference)
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = 'name: nf-core/fastqc\nversion: 1.0.0\n'
        tempDir.resolve(RemoteModuleResolver.PIPELINE_INFO_FILE).text = 'name=a/b\n'

        when:
        def result = new ModuleResolver(tempDir).resolve(reference, null, false)

        then:
        result == moduleDir.resolve('main.nf')
        !appender.list*.formattedMessage.any { it.contains('no registry origin') }

        cleanup:
        logger.detachAppender(appender)
        logger.setLevel(level)
    }
}
