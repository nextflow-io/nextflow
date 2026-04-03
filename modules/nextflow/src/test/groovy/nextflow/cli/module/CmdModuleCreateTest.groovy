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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdModuleCreateTest extends Specification {

    @TempDir
    Path tempDir

    def 'should generate main.nf content'() {
        when:
        def content = CmdModuleCreate.mainNf('myorg', 'hello')

        then:
        content.contains("Module: myorg/hello")
        content.contains("process HELLO")
    }

    def 'should generate README.md content'() {
        when:
        def content = CmdModuleCreate.readmeMd('myorg', 'hello')

        then:
        content.contains('# myorg/hello')
        content.contains('## Summary')
        content.contains('## Get started')
        content.contains('## Input / Output / Dependencies')
        content.contains('## License')
        content.contains("include { HELLO } from 'myorg/hello'")
    }

    def 'should reject invalid args'() {
        given:
        def cmd = new CmdModuleCreate(args: ['foo', 'bar'])

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should reject invalid format'() {
        given:
        def cmd = new CmdModuleCreate(args: ['invalid'])

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail if directory already exists'() {
        given:
        def moduleDir = tempDir.resolve('modules/myorg/hello')
        Files.createDirectories(moduleDir)
        and:
        def cmd = Spy(CmdModuleCreate) {
            modulesBase() >> tempDir.resolve('modules')
        }

        when:
        cmd.createModule('myorg', 'hello')

        then:
        thrown(AbortOperationException)
    }

    def 'should create module directory and files'() {
        given:
        def cmd = Spy(CmdModuleCreate) {
            modulesBase() >> tempDir.resolve('modules')
        }

        when:
        cmd.createModule('testorg', 'testmod')

        then:
        def moduleDir = tempDir.resolve('modules/testorg/testmod')
        Files.exists(moduleDir.resolve('main.nf'))
        Files.exists(moduleDir.resolve('README.md'))
        Files.exists(moduleDir.resolve('meta.yml'))
        Files.exists(moduleDir.resolve('.module-info'))
        moduleDir.resolve('main.nf').text.contains('process TESTMOD')
        moduleDir.resolve('README.md').text.contains('# testorg/testmod')
        moduleDir.resolve('meta.yml').text.contains('name: testorg/testmod')
        moduleDir.resolve('meta.yml').text.contains('version: 1.0.0')
    }

    def 'should create .module-info file'() {
        given:
        def cmd = Spy(CmdModuleCreate) {
            modulesBase() >> tempDir.resolve('modules')
        }

        when:
        cmd.createModule('myorg', 'hello')

        then:
        def moduleDir = tempDir.resolve('modules/myorg/hello')
        Files.exists(moduleDir.resolve('.module-info'))
    }
}
