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

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import nextflow.util.CacheHelper
import spock.lang.Specification

class SpecTaskHasherExplainTest extends Specification {

    private HashContext ctx(String source, String conda) {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) {
                getHashMode() >> CacheHelper.HashMode.STANDARD
            }
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
        }
        def task = Mock(TaskRun) {
            getSource() >> source
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [:]
            getCondaEnv() >> Paths.get(conda)
            getSpackEnv() >> null
            hasStubBlock() >> false
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [:]
        helper.getTaskBinEntries(_) >> []
        return new HashContext(task, helper)
    }

    def 'explain reports one digest per emitting key, in spec order'() {
        when:
        def explained = new SpecTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4).explain()

        then: 'keys that emitted nothing are absent'
        explained.keySet() as List == [
            HashKey.SESSION_ID, HashKey.PROCESS_NAME, HashKey.TASK_SOURCE, HashKey.CONDA
        ]
    }

    def 'only the differing key changes digest between two tasks'() {
        given:
        def a = new SpecTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4).explain()
        def b = new SpecTaskHasher(ctx('echo CHANGED', 'env.yml'), StdSpecs.STD_V4).explain()

        expect:
        a[HashKey.TASK_SOURCE] != b[HashKey.TASK_SOURCE]
        and:
        (a.keySet() - HashKey.TASK_SOURCE).every { a[it] == b[it] }
    }

    def 'explain does not change the computed hash'() {
        given:
        def hasher = new SpecTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4)
        def before = hasher.compute()

        when:
        hasher.explain()

        then:
        hasher.compute() == before
    }
}
