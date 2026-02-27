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
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for ModulePublish command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModulePublishTest extends Specification {

    @TempDir
    Path tempDir

    def 'should validate module structure' () {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)

        // Create required files
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('README.md').text = '# Test Module'
        moduleDir.resolve('meta.yml').text = '''
name: test/module
version: 1.0.0
description: Test module
license: MIT
'''

        and:
        def cmd = new ModulePublish()
        cmd.dryRun = true
        cmd.args = [moduleDir.toString()]

        when:
        def errors = cmd.invokeMethod('validateModuleStructure', moduleDir)

        then:
        errors.isEmpty()
    }

    def 'should detect missing required files' () {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)

        // Only create main.nf, missing meta.yaml and README.md
        moduleDir.resolve('main.nf').text = 'process TEST { }'

        and:
        def cmd = new ModulePublish()

        when:
        def errors = cmd.invokeMethod('validateModuleStructure', moduleDir)

        then:
        errors.size() == 2
        errors.any { it.contains('meta.yml') }
        errors.any { it.contains('README.md') }
    }

    def 'should detect oversized module' () {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)

        // Create required files
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('README.md').text = '# Test Module'
        moduleDir.resolve('meta.yml').text = '''
name: test/module
version: 1.0.0
description: Test module
license: MIT
'''

        // Create a large file (>1MB)
        def largeFile = moduleDir.resolve('large-file.txt')
        def content = 'x' * (1024 * 1024 + 1000) // 1MB + 1000 bytes
        Files.writeString(largeFile, content)

        and:
        def cmd = new ModulePublish()

        when:
        def errors = cmd.invokeMethod('validateModuleStructure', moduleDir)

        then:
        errors.size() == 1
        errors[0].contains('1MB limit')
    }

    def 'should succeed in dry-run mode without authentication' () {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)

        // Create required files
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('README.md').text = '# Test Module'
        moduleDir.resolve('meta.yml').text = '''
name: test/module
version: 1.0.0
description: Test module
license: MIT
'''

        and:
        def cmd = new ModulePublish()
        def launcher = new Launcher()
        launcher.options = [:]
        cmd.launcher = launcher
        cmd.args = [moduleDir.toString()]
        cmd.dryRun = true

        when:
        cmd.run()

        then:
        noExceptionThrown()
    }
}
