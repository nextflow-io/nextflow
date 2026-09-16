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
package nextflow.processor.hash

import groovy.json.JsonSlurper
import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import nextflow.util.CacheHelper
import spock.lang.Specification

class BaseTaskHasherDumpTest extends Specification {

    def 'dumpJson names every entry and identifies the spec'() {
        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
            getStubBlock() >> null
            getHashMode() >> CacheHelper.HashMode.STANDARD
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo a'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [:]
            getCondaEnv() >> null
            getSpackEnv() >> null
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [:]
        helper.getTaskBinEntries(_) >> []

        when:
        def json = new JsonSlurper().parseText(
            new BaseTaskHasher(new HashContext(task, helper), StdSpecs.STD_V4).dumpJson()) as List

        then:
        json[0].spec == 'std/v4'
        json[0].fingerprint == StdSpecs.STD_V4.fingerprint()
        and:
        json[1..-1]*.key == ['SESSION_ID', 'PROCESS_NAME', 'TASK_SOURCE']
        and:
        json[1..-1].every { it.hash instanceof String && it.hash.length() > 0 }
    }

    def 'dumpLegacy produces non-empty output with header and entries'() {
        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
            getStubBlock() >> null
            getHashMode() >> CacheHelper.HashMode.STANDARD
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo a'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [:]
            getCondaEnv() >> null
            getSpackEnv() >> null
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [:]
        helper.getTaskBinEntries(_) >> []

        when:
        def output = new BaseTaskHasher(new HashContext(task, helper), StdSpecs.STD_V4).dumpLegacy()

        then:
        output.size() > 0
        and:
        output.contains('PIPE:FOO')
        output.contains('cache hash:')
        output.contains('mode:')
        output.contains('entries:')
        and:
        output.split('\n').findAll { it.startsWith('  ') }.size() > 0
    }
}
