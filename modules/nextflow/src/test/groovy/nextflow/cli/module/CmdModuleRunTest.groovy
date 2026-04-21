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

import nextflow.cli.CliOptions
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for CmdModuleRun source resolution
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CmdModuleRunTest extends Specification {

    @TempDir
    Path tempDir

    // --- Source classification (data-driven) ---

    def 'isLocalModule classifies "#source" as local=#expected'() {
        expect:
        new CmdModuleRun().isLocalModule(source) == expected

        where:
        source                  || expected
        '/absolute/path'        || true
        './relative/path'       || true
        '../parent/path'        || true
        '/tmp/module/main.nf'   || true
        'nf-core/test-module'   || false
        'scope/name'            || false
        'single-word'           || false
    }

    // --- Local resolution semantics ---

    def 'local module directory resolves to main.nf'() {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = '// local module'

        when:
        def result = new CmdModuleRun().resolveLocalModule(moduleDir.toString())

        then:
        result == moduleDir.resolve('main.nf')
    }

    def 'local module file resolves directly'() {
        given:
        def scriptFile = tempDir.resolve('custom.nf')
        scriptFile.text = '// custom script'

        when:
        def result = new CmdModuleRun().resolveLocalModule(scriptFile.toString())

        then:
        result == scriptFile
    }

    // --- Source dispatch interactions ---

    def 'local module directory dispatches to local resolver and bypasses remote resolver'() {
        given:
        def moduleDir = tempDir.resolve('my-module')
        Files.createDirectories(moduleDir)
        moduleDir.resolve('main.nf').text = '// local module'

        and:
        def resolved = moduleDir.resolve('main.nf')
        def cmd = Spy(CmdModuleRun)
        cmd.args = [moduleDir.toString()]

        when:
        def result = cmd.resolveModuleSource()

        then:
        result == resolved
        1 * cmd.resolveLocalModule(moduleDir.toString()) >> resolved
        0 * cmd.resolveRemoteModule(_, _)
    }

    def 'local module file dispatches to local resolver and bypasses remote resolver'() {
        given:
        def scriptFile = tempDir.resolve('custom.nf')
        scriptFile.text = '// custom script'

        and:
        def cmd = Spy(CmdModuleRun)
        cmd.args = [scriptFile.toString()]

        when:
        def result = cmd.resolveModuleSource()

        then:
        result == scriptFile
        1 * cmd.resolveLocalModule(scriptFile.toString()) >> scriptFile
        0 * cmd.resolveRemoteModule(_, _)
    }

    def 'registry source dispatches to remote resolver and bypasses local resolver'() {
        given:
        def resolved = tempDir.resolve('modules/nf-core/test-module/main.nf')
        def cmd = Spy(CmdModuleRun)
        cmd.args = ['nf-core/test-module']
        cmd.launcher = Mock(Launcher) { getOptions() >> new CliOptions() }

        when:
        def result = cmd.resolveModuleSource()

        then:
        result == resolved
        1 * cmd.resolveRemoteModule('nf-core/test-module', null) >> resolved
        0 * cmd.resolveLocalModule(_)
    }

    // --- Error cases ---

    def 'no arguments throws AbortOperationException'() {
        given:
        def cmd = new CmdModuleRun()
        cmd.args = []

        when:
        cmd.resolveModuleSource()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Module name/path not provided')
    }

    def 'non-existent local path throws AbortOperationException'() {
        given:
        def cmd = new CmdModuleRun()
        cmd.args = [tempDir.resolve('does-not-exist').toString()]
        cmd.root = tempDir

        when:
        cmd.resolveModuleSource()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Invalid module path')
    }

    def 'invalid module reference format throws AbortOperationException'() {
        given:
        def cmd = new CmdModuleRun()
        cmd.args = ['invalid-module']
        cmd.root = tempDir
        cmd.launcher = Mock(Launcher) { getOptions() >> new CliOptions() }

        when:
        cmd.resolveModuleSource()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Invalid module reference')
    }
}
