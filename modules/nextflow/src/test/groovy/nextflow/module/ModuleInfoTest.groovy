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

import spock.lang.Specification

/**
 * Test suite for ModuleInfo
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleInfoTest extends Specification {

    Path tempDir

    def setup() {
        tempDir = Files.createTempDirectory('nf-module-info-test-')
    }

    def cleanup() {
        tempDir?.deleteDir()
    }

    def 'should save and load a single property'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        ModuleInfo.save(moduleDir, 'checksum', 'abc123')

        then:
        ModuleInfo.load(moduleDir, 'checksum') == 'abc123'
    }

    def 'should create .module-info file on save'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        ModuleInfo.save(moduleDir, 'version', '1.0.0')

        then:
        Files.exists(moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE))
    }

    def 'should update existing property without affecting others'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        ModuleInfo.save(moduleDir, 'checksum', 'original')
        ModuleInfo.save(moduleDir, 'registryUrl', 'http://registry.com')

        when:
        ModuleInfo.save(moduleDir, 'checksum', 'updated')

        then:
        ModuleInfo.load(moduleDir, 'checksum') == 'updated'
        ModuleInfo.load(moduleDir, 'registryUrl') == 'http://registry.com'
    }

    def 'should return null when loading property from non-existent file'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        def result = ModuleInfo.load(moduleDir, 'checksum')

        then:
        result == null
    }

    def 'should return null when loading non-existent property'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        ModuleInfo.save(moduleDir, 'registryUrl', 'http://registry.com')

        when:
        def result = ModuleInfo.load(moduleDir, 'missing-property')

        then:
        result == null
    }

    def 'should save and load multiple properties via map'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        def props = [checksum: 'abc123', registryUrl: 'http://registry.com', author: 'test']

        when:
        ModuleInfo.save(moduleDir, props)

        then:
        ModuleInfo.load(moduleDir, 'checksum') == 'abc123'
        ModuleInfo.load(moduleDir, 'registryUrl') == 'http://registry.com'
        ModuleInfo.load(moduleDir, 'author') == 'test'
    }

    def 'should merge map properties with existing ones'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        ModuleInfo.save(moduleDir, 'existing', 'value')

        when:
        ModuleInfo.save(moduleDir, [newprop: 'newval'])

        then:
        ModuleInfo.load(moduleDir, 'existing') == 'value'
        ModuleInfo.load(moduleDir, 'newprop') == 'newval'
    }

    def 'should do nothing when saving null map'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        ModuleInfo.save(moduleDir, (Map) null)

        then:
        !Files.exists(moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE))
    }

    def 'should do nothing when saving empty map'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        ModuleInfo.save(moduleDir, [:])

        then:
        !Files.exists(moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE))
    }

    def 'should load all properties as map'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        ModuleInfo.save(moduleDir, [checksum: 'abc123', registryUrl: 'http://registry.com'])

        when:
        def result = ModuleInfo.load(moduleDir)

        then:
        result['checksum'] == 'abc123'
        result['registryUrl'] == 'http://registry.com'
    }

    def 'should return empty map when loading all from non-existent file'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        def result = ModuleInfo.load(moduleDir)

        then:
        result == [:]
    }
}
