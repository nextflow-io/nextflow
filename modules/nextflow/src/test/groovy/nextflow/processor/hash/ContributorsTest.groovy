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
import nextflow.processor.Architecture
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

class ContributorsTest extends Specification {

    private HashContext ctxFor(TaskRun task, TaskHasher helper) {
        return new HashContext(task, helper)
    }

    def 'absent optional keys emit nothing, not null'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            isContainerEnabled() >> false
            getOutputEvals() >> [:]
            getCondaEnv() >> null
            getSpackEnv() >> null
            getConfig() >> config
        }
        def helper = Mock(TaskHasher)
        def ctx = ctxFor(task, helper)

        expect:
        Contributors.CONTAINER.emit(ctx) == []
        Contributors.EVAL_OUTPUTS_RAW_MAP.emit(ctx) == []
        Contributors.CONDA.emit(ctx) == []
        Contributors.SPACK_AND_ARCH.emit(ctx) == []
        Contributors.ENV_MODULES.emit(ctx) == []
    }

    def 'spack emits arch only when spack itself is set'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def config = Mock(TaskConfig) { getArchitecture() >> new Architecture('linux/amd64') }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getSpackEnv() >> spack
            getConfig() >> config
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.SPACK_AND_ARCH.emit(ctx) == expected

        where:
        spack                       | expected
        null                        | []
        // getSpackEnv() returns Path and getArchitecture() returns Architecture, not
        // String — a plain String stub here makes Spock's mock response coercion throw
        // GroovyCastException trying to cast it to the mocked method's declared type.
        Paths.get('env.yaml')       | [Paths.get('env.yaml'), new Architecture('linux/amd64')]
    }

    def 'inputs emit a name and a value per input, in order'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getInputs() >> [
                (Mock(nextflow.script.params.InParam) { getName() >> 'x' }): 1,
                (Mock(nextflow.script.params.InParam) { getName() >> 'y' }): 2
            ]
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.INPUTS_RAW.emit(ctx) == ['x', 1, 'y', 2]
    }

    def 'the two eval forms differ, and the derived form is the pre-7575 string'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def evals = [beta: 'echo b', alpha: 'echo a']
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getOutputEvals() >> evals
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.EVAL_OUTPUTS_RAW_MAP.emit(ctx) == ['eval_outputs', evals]
        and:
        Contributors.EVAL_OUTPUTS_DERIVED_STRING.emit(ctx) == ['eval_outputs', 'alpha=echo a\nbeta=echo b']
    }

    def 'canonical names are distinct across all contributors'() {
        given:
        def all = [
            Contributors.SESSION_ID, Contributors.PROCESS_NAME, Contributors.TASK_SOURCE,
            Contributors.CONTAINER, Contributors.INPUTS_RAW, Contributors.EVAL_OUTPUTS_RAW_MAP,
            Contributors.EVAL_OUTPUTS_DERIVED_STRING, Contributors.SCRIPT_VARS,
            Contributors.BIN_ENTRIES, Contributors.MODULE_BUNDLE, Contributors.ENV_MODULES,
            Contributors.CONDA, Contributors.SPACK_AND_ARCH, Contributors.STUB_MARKER
        ]

        expect:
        all*.canonicalName().toSet().size() == all.size()
    }
}
