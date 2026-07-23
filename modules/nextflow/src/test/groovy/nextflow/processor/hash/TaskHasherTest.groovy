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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Session
import nextflow.util.CacheHelper
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import nextflow.script.bundle.ResourcesBundle
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHasherTest extends Specification {

    def 'should compute unique task hash' () {

        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            getBinEntries() >> [:]
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
        }
        and:
        def task = Mock(TaskRun) {
            getSource() >> 'hello world'
            isContainerEnabled() >> false
            getContainer() >> null
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
        }
        and:
        def hasher = Spy(new TaskHasherV1(task))

        when:
        def uuid1 = hasher.compute()
        def uuid2 = hasher.compute()
        then:
        hasher.getTaskGlobalVars() >>> [
            [foo:'a', bar:'b'],
            [bar:'b', foo:'a']
        ]
        and: 'global vars should not affect task hash'
        uuid1 == uuid2
    }

    def 'should include referenced bin files in the task hash' () {

        given:
        def session = Mock(Session) {
            getBinEntries() >> [
                'foo.sh': Path.of('/some/path/foo.sh'),
                'bar.sh': Path.of('/some/path/bar.sh')
            ]
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
        }
        def task = new TaskRun()
        task.processor = processor

        when:
        def result = task.getTaskBinEntries('var=x foo.sh')
        then:
        result.size() == 1
        result.contains(Path.of('/some/path/foo.sh'))

        when:
        result = task.getTaskBinEntries('echo $(foo.sh); bar.sh')
        then:
        result.size() == 2
        result.contains(Path.of('/some/path/foo.sh'))
        result.contains(Path.of('/some/path/bar.sh'))
    }

    def 'should include module bundle fingerprint in the task hash' () {

        given: 'two module bundles with the same fingerprint, then one with a different fingerprint'
        def bundleA1 = Mock(ResourcesBundle) { asBoolean() >> true; hasEntries() >> true; fingerprint() >> 'aaaaaaaa' }
        def bundleA2 = Mock(ResourcesBundle) { asBoolean() >> true; hasEntries() >> true; fingerprint() >> 'aaaaaaaa' }
        def bundleB = Mock(ResourcesBundle) { asBoolean() >> true; hasEntries() >> true; fingerprint() >> 'bbbbbbbb' }
        and: 'a task whose process resolves those bundles across successive hash computations'
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            getBinEntries() >> [:]
            enableModuleBinaries() >> true
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >>> [bundleA1, bundleA2, bundleB]
        }
        def task = Mock(TaskRun) {
            getSource() >> 'hello world'
            isContainerEnabled() >> false
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
        }
        def hasher = Spy(new TaskHasherV2(task))
        hasher.getTaskGlobalVars() >> [:]

        when:
        def hashA1 = hasher.compute()
        def hashA2 = hasher.compute()
        def hashB = hasher.compute()

        then: 'the same bundle fingerprint yields the same hash'
        hashA1 == hashA2
        and: 'a different bundle fingerprint invalidates the cache'
        hashA1 != hashB
    }

    def 'should hash bin entries only in V3' () {

        given: 'a referenced bin script backed by a real file'
        def tmp = Files.createTempFile('foo', '.sh')
        tmp.text = 'echo A'
        and:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) { getHashMode() >> CacheHelper.HashMode.STANDARD }
        }
        def task = Mock(TaskRun) {
            getSource() >> 'foo.sh'
            isContainerEnabled() >> false
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
            getTaskBinEntries('foo.sh') >> [tmp]
        }
        and:
        def v2 = Spy(new TaskHasherV2(task)); v2.getTaskGlobalVars() >> [:]
        def v3 = Spy(new TaskHasherV3(task)); v3.getTaskGlobalVars() >> [:]

        when: 'the referenced bin script content changes between computations'
        def v2before = v2.compute()
        def v3before = v3.compute()
        tmp.text = 'echo B'
        def v2after = v2.compute()
        def v3after = v3.compute()

        then: 'V3 tracks the bin script content'
        v3before != v3after
        and: 'V2 ignores bin entries'
        v2before == v2after

        cleanup:
        Files.deleteIfExists(tmp)
    }

    def 'should compute hash entries for eval outputs'() {

        when:
        def result1 = AbstractTaskHasher.computeEvalOutputCommands([
            'nxf_out_eval_2': 'echo "value2"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_3': 'echo "value3"'
        ])

        def result2 = AbstractTaskHasher.computeEvalOutputCommands([
            'nxf_out_eval_3': 'echo "value3"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_2': 'echo "value2"'
        ])

        then:
        result1 == result2
        result1 == 'nxf_out_eval_1=echo "value1"\nnxf_out_eval_2=echo "value2"\nnxf_out_eval_3=echo "value3"'
    }
}
