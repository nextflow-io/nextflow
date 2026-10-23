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

    // a plugin compiled before `exec(String,String,List,Map)` was added
    static class LegacyCommand implements PluginExecAware {
        List<String> received
        int exec(String pluginId, String cmd, List<String> args) { received = [pluginId, cmd]; return 0 }
    }

    def cleanup() {
        SysEnv.pop()
    }

    def 'should create the session with the given config' () {
        given:
        SysEnv.push([:])
        def command = new TestCommand()

        when:
        command.exec('nf-test', 'test', [], [workDir: '/my/work', azure: [storage: [accountName: 'my-account']]])
        then:
        command.config.workDir == '/my/work'
        command.config.azure.storage.accountName == 'my-account'
    }

    def 'should create the session with the env config' () {
        given:
        SysEnv.push([NXF_WORK: '/env/work'])
        def command = new TestCommand()

        when:
        command.exec('nf-test', 'test', [])
        then:
        command.config.workDir == '/env/work'
    }

    def 'should fallback to legacy exec method' () {
        given:
        SysEnv.push([:])
        def command = new LegacyCommand()

        when:
        def ret = command.exec('nf-test', 'test', [], [workDir: '/my/work'])
        then:
        ret == 0
        command.received == ['nf-test', 'test']
    }

}
