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

package nextflow.cli

import java.nio.file.Files

import nextflow.SysEnv
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginAbstractExecTest extends Specification {

    static class TestCommand implements PluginAbstractExec {
        Map config
        List<String> getCommands() { ['test'] }
        int exec(String cmd, List<String> args) { config = getSession().config; return 0 }
    }

    def cleanup() {
        SysEnv.pop()
    }

    def 'should create the session with the config file and the work dir env override' () {
        given:
        def folder = Files.createTempDirectory('test')
        def configFile = folder.resolve('nextflow.config')
        configFile.text = '''
            workDir = "$WORK_DIR"
            azure.storage.accountName = 'my-account'
            '''.replace('$WORK_DIR', folder.resolve('config-work').toString())
        and:
        SysEnv.push(ENV + [NXF_CONFIG_FILE: configFile.toString()])
        def command = new TestCommand()

        when:
        command.exec('nf-test', 'test', [])
        then:
        command.config.workDir == (WORK ?: folder.resolve('config-work').toString())
        command.config.azure.storage.accountName == 'my-account'

        cleanup:
        folder?.deleteDir()

        where:
        ENV                     | WORK
        [:]                     | null
        [NXF_WORK: '/env/work'] | '/env/work'
    }

    def 'should not resolve config includes' () {
        given:
        def folder = Files.createTempDirectory('test')
        def configFile = folder.resolve('nextflow.config')
        configFile.text = '''
            azure.storage.accountName = 'my-account'
            includeConfig 'https://localhost:1/single-use.config'
            '''
        and:
        SysEnv.push([NXF_CONFIG_FILE: configFile.toString()])
        def command = new TestCommand()

        when:
        command.exec('nf-test', 'test', [])
        then:
        command.config.azure.storage.accountName == 'my-account'

        cleanup:
        folder?.deleteDir()
    }

}
