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
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

class HashContextTest extends Specification {

    def 'exposes task, processor and session, and delegates the helpers'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) {
            getSession() >> session
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getSource() >> 'echo hello'
        }
        and:
        def helper = Mock(TaskHasher) {
            getTaskGlobalVars() >> [foo: 'a']
            getTaskBinEntries('echo hello') >> [Paths.get('/bin/x.sh')]
        }
        and:
        def ctx = new HashContext(task, helper)

        expect:
        ctx.task.is(task)
        ctx.processor.is(processor)
        ctx.session.is(session)
        ctx.globalVars() == [foo: 'a']
        ctx.binEntries() == [Paths.get('/bin/x.sh')]
    }

    def 'the vocabulary has exactly the keys master hashes'() {
        expect:
        HashKey.values()*.name() as Set == [
            'SESSION_ID', 'PROCESS_NAME', 'TASK_SOURCE', 'CONTAINER', 'INPUTS',
            'EVAL_OUTPUTS', 'SCRIPT_VARS', 'BIN_ENTRIES', 'MODULE_BUNDLE',
            'ENV_MODULES', 'CONDA', 'SPACK', 'STUB_MARKER'
        ] as Set
    }
}
