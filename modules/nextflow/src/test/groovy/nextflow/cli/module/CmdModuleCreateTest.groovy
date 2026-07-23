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

    def 'should generate a workflow module main.nf with -kind Workflow'() {
        when:
        def content = CmdModuleCreate.mainNf('myorg', 'hello', 'Workflow')

        then:
        content.contains("Workflow module: myorg/hello")
        content.contains("workflow HELLO")
        content.contains("take:")
        content.contains("emit:")
        !content.contains("process ")
    }

    def 'should generate a workflow module meta.yml with kind Workflow'() {
        when:
        def content = CmdModuleCreate.metaYml('myorg', 'hello', 'Workflow')

        then:
        content.contains("name: myorg/hello")
        content.contains("kind: Workflow")

        and: 'the untyped scaffold take/emit are documented as channels'
        content.contains("input:")
        content.contains("name: ch_input")
        content.contains("type: channel")
        content.contains("name: output")

        and: 'untyped workflows require the baseline Nextflow version'
        content.contains('nextflow: ">=24.04.0"')
    }

    def 'generated untyped workflow scaffold passes validation'() {
        given:
        def moduleDir = tempDir.resolve('modules/myorg/hello')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = CmdModuleCreate.mainNf('myorg', 'hello', 'Workflow')
        moduleDir.resolve('meta.yml').text = CmdModuleCreate.metaYml('myorg', 'hello', 'Workflow')
        moduleDir.resolve('README.md').text = '# hello\n'
        and: 'a permissive schema'
        def schema = tempDir.resolve('schema.json')
        schema.text = '{"$schema":"https://json-schema.org/draft/2020-12/schema","type":"object",' +
            '"properties":{"name":{"type":"string"},"description":{"type":"string"}},"required":["name","description"]}'

        expect: 'no errors -- interface counts match take/emit and there are no TODO placeholders'
        nextflow.module.ModuleValidator.validate(moduleDir, schema.toString()).isEmpty()
    }

    def 'should default to a process module'() {
        expect:
        CmdModuleCreate.mainNf('myorg', 'hello').contains("process HELLO")
        !CmdModuleCreate.metaYml('myorg', 'hello').contains("kind:")
    }

    def 'should generate a typed process module with -typed'() {
        when:
        def content = CmdModuleCreate.mainNf('myorg', 'hello', 'Process', true)

        then:
        content.contains("nextflow.enable.types = true")
        content.contains("process HELLO")
        content.contains("greeting: String")               // typed input
        content.contains("message: String = stdout()")     // named typed output bound to stdout
        content.contains("script:")                          // keeps a shell script block
        !content.contains("val greeting")                    // not the untyped form
        !content.contains("exec:")

        and:
        def meta = CmdModuleCreate.metaYml('myorg', 'hello', 'Process', true)
        meta.contains("name: message")
        meta.contains('nextflow: ">=25.10.0"')   // typed processes require Nextflow 25.10.0+
    }

    def 'generated typed workflow scaffold parses as a workflow'() {
        given:
        def nf = tempDir.resolve('main.nf')
        nf.text = CmdModuleCreate.mainNf('myorg', 'hello', 'Workflow', true)

        expect:
        nextflow.module.ModuleSpecFactory.definesWorkflow(nf)
    }

    def 'generated typed process scaffold parses cleanly'() {
        given:
        def nf = tempDir.resolve('main.nf')
        nf.text = CmdModuleCreate.mainNf('myorg', 'hello', 'Process', true)

        expect:
        // parses without error (no workflow defined -> false, but must not throw)
        !nextflow.module.ModuleSpecFactory.definesWorkflow(nf)
    }

    def 'should generate a typed workflow module with -kind Workflow -typed'() {
        when:
        def content = CmdModuleCreate.mainNf('myorg', 'hello', 'Workflow', true)

        then:
        content.contains("nextflow.enable.types = true")
        content.contains("workflow HELLO")
        content.contains("greeting: String")            // typed take
        content.contains("result: String = message")    // typed emit
        !content.contains("process ")

        and: 'the meta.yml derives input/output from the take/emit'
        def meta = CmdModuleCreate.metaYml('myorg', 'hello', 'Workflow', true)
        meta.contains("kind: Workflow")
        meta.contains("input:")
        meta.contains("name: greeting")
        meta.contains("output:")
        meta.contains("name: result")

        and: 'typed workflows require Nextflow 26.04.0+'
        meta.contains('nextflow: ">=26.04.0"')
    }

    def 'should translate special chars in process name to underscore and uppercase'() {
        expect:
        CmdModuleCreate.mainNf('myorg', name).contains("process ${expected}")

        where:
        name            | expected
        'hello'         | 'HELLO'
        'hello-world'   | 'HELLO_WORLD'
        'hello.world'   | 'HELLO_WORLD'
        'hello world'   | 'HELLO_WORLD'
        'my-tool_v2'    | 'MY_TOOL_V2'
        'foo/bar'       | 'FOO_BAR'
    }

    def 'should reject invalid namespace'() {
        given:
        def cmd = new CmdModuleCreate(args: ["${namespace}/hello".toString()])

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('namespace')

        where:
        namespace << ['.hidden', '..', 'has space', '../evil', ]
    }

    def 'should reject invalid name'() {
        given:
        def cmd = new CmdModuleCreate(args: ["myorg/${name}".toString()])

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('name')

        where:
        name << ['.hidden', '..', 'has space', 'sub/.hidden', 'sub/..', 'sub/has space']
    }

    def 'should accept valid namespace and name'() {
        given:
        def cmd = Spy(CmdModuleCreate, constructorArgs: [[args: ["${namespace}/${name}".toString()]]]) {
            modulesBase() >> tempDir.resolve('modules')
        }

        when:
        cmd.run()

        then:
        noExceptionThrown()
        Files.exists(tempDir.resolve("modules/${namespace}/${name}/main.nf"))

        where:
        namespace   | name
        'myorg'     | 'hello'
        'my-org'    | 'my-tool'
        'org123'    | 'tool.v2'
        'My_Org'    | 'My_Tool'
        'myorg'     | 'submodule/hello'
        'myorg'     | 'sub/nested/tool'
    }

    def 'should generate README.md content'() {
        when:
        def content = CmdModuleCreate.readmeMd('myorg', 'hello')

        then:
        content.contains('# myorg/hello')
        content.contains('## Summary')
        content.contains('## Get started')
        content.contains('## Dependencies')
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
