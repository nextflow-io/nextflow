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
import nextflow.module.ModuleValidator
import spock.lang.Specification
import spock.lang.TempDir

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdModuleValidateTest extends Specification {

    @TempDir
    Path tempDir

    private Path createValidModule(String namespace='myorg', String name='hello') {
        def moduleDir = tempDir.resolve("modules/$namespace/$name")
        Files.createDirectories(moduleDir)

        moduleDir.resolve('main.nf').text = '''\
            process HELLO {
                input:
                val greeting

                output:
                stdout

                script:
                """
                echo '${greeting}'
                """
            }
            '''.stripIndent()

        moduleDir.resolve('meta.yml').text = """\
            name: ${namespace}/${name}
            version: 1.0.0
            description: A test module
            license: Apache-2.0
            input:
              - name: greeting
                type: string
                description: A greeting string
            """.stripIndent()

        moduleDir.resolve('README.md').text = "# ${namespace}/${name}\n"
        moduleDir.resolve('.module-info').text = ''

        return moduleDir
    }

    def 'should pass validation for valid module'() {
        given:
        def moduleDir = createValidModule()

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.isEmpty()
    }

    def 'should fail when main.nf is missing'() {
        given:
        def moduleDir = createValidModule()
        Files.delete(moduleDir.resolve('main.nf'))

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.any { it.contains('main.nf') }
    }

    def 'should fail when meta.yml is missing'() {
        given:
        def moduleDir = createValidModule()
        Files.delete(moduleDir.resolve('meta.yml'))

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.any { it.contains('meta.yml') }
    }

    def 'should fail when README.md is missing'() {
        given:
        def moduleDir = createValidModule()
        Files.delete(moduleDir.resolve('README.md'))

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.any { it.contains('README.md') }
    }

    def 'should fail when meta.yml has missing required fields'() {
        given:
        def moduleDir = createValidModule()
        moduleDir.resolve('meta.yml').text = '''\
            name: myorg/hello
            version: 1.0.0
            '''.stripIndent()

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.any { it.contains('description') }
        errors.any { it.contains('license') }
    }

    def 'should fail when meta.yml has invalid version format'() {
        given:
        def moduleDir = createValidModule()
        moduleDir.resolve('meta.yml').text = '''\
            name: myorg/hello
            version: bad
            description: test
            license: MIT
            '''.stripIndent()

        when:
        def errors = ModuleValidator.validate(moduleDir)

        then:
        errors.any { it.contains('version') }
    }

    def 'should fail for non-existent module directory'() {
        when:
        def errors = ModuleValidator.validateStructure(tempDir.resolve('nonexistent'))

        then:
        errors.any { it.contains('does not exist') }
    }

    def 'should throw for invalid argument'() {
        given:
        def cmd = new CmdModuleValidate(args: ['not-a-module'], root: tempDir)

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }
}
