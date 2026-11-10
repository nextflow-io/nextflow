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

package nextflow.processor

import java.nio.file.Path

import nextflow.Session
import nextflow.file.FileHolder
import nextflow.script.ProcessConfig
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.v2.ProcessInput
import nextflow.util.CacheHelper
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
        def hasher = Spy(new TaskHasher(task))

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
        def task = Mock(TaskRun) {
            getProcessor() >> processor
        }
        def hasher = new TaskHasher(task)

        when:
        def result = hasher.getTaskBinEntries('var=x foo.sh')
        then:
        result.size() == 1
        result.contains(Path.of('/some/path/foo.sh'))

        when:
        result = hasher.getTaskBinEntries('echo $(foo.sh); bar.sh')
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
        def hasher = Spy(new TaskHasher(task))
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

    def 'should get task directive vars' () {
        given:
        def processor = Spy(TaskProcessor) {
            getConfig() >> Mock(ProcessConfig)
        }
        and:
        def config = new TaskConfig()
        config.cpus = 4
        config.ext.alpha = 'AAAA'
        config.ext.delta = { foo }
        config.ext.omega = "${-> bar}"
        config.context = [foo: 'DDDD', bar: 'OOOO']
        and:
        def task = Mock(TaskRun) {
            getConfig() >> config
            getProcessor() >> processor
            getVariableNames() >> { [ 'task.cpus', 'task.ext.alpha', 'task.ext.delta', 'task.ext.omega' ] as Set }
        }

        when:
        def result = new TaskHasher(task).getTaskExtensionDirectiveVars()
        then:
        result == [
            'task.ext.alpha': 'AAAA',
            'task.ext.delta': 'DDDD',
            'task.ext.omega': 'OOOO',
        ]
    }

    def 'should hash eval outputs independently of map order'() {
        given:
        // same entries, different insertion order and different map implementation
        def evals1 = new LinkedHashMap(['nxf_out_eval_0': 'echo one', 'nxf_out_eval_1': 'echo two'])
        def evals2 = new HashMap(['nxf_out_eval_1': 'echo two', 'nxf_out_eval_0': 'echo one'])

        expect:
        CacheHelper.hasher(evals1).hash() == CacheHelper.hasher(evals2).hash()
    }

    def 'should hash eval outputs depending on the name-command binding'() {
        given:
        // same names and same commands, but paired the other way round
        def evals1 = ['nxf_out_eval_0': 'echo one', 'nxf_out_eval_1': 'echo two']
        def evals2 = ['nxf_out_eval_0': 'echo two', 'nxf_out_eval_1': 'echo one']

        expect:
        CacheHelper.hasher(evals1).hash() != CacheHelper.hasher(evals2).hash()
    }

    private List<Object> stageKeys(Map params) {
        final processor = Mock(TaskProcessor) {
            getSession() >> Mock(Session)
        }
        final task = Mock(TaskRun) {
            getProcessor() >> processor
            getInputEnv() >> params.env
            getStdin() >> params.stdin
            getStageAsFiles() >> params.staged
        }
        final keys = new ArrayList<Object>()
        new TaskHasher(task).addStageDirectives(keys)
        return keys
    }

    private String stageHash(Map params, CacheHelper.HashMode mode = CacheHelper.HashMode.STANDARD) {
        return CacheHelper.hasher(stageKeys(params), mode).hash().toString()
    }

    def 'should not add stage keys when the stage block is not used' () {
        expect:
        stageKeys([:]) == []
    }

    def 'should hash stage directives to a stable value' () {
        given:
        def staged = new TaskStageAsFiles('query.fa', [new FileHolder(Path.of('/data/a.fa')).withName('query.fa')], 'fasta')

        expect:
        stageHash(env: [FOO: 'one'], stdin: 'two', staged: [staged]) == '5d722237e019c1bd176c8cf169bf1ece'
    }

    def 'should keep a stable hash for typed tasks without stage directives' () {
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
        def task = Mock(TaskRun) {
            getSource() >> 'echo hello'
            isContainerEnabled() >> false
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
            hasTypedInputsOutputs() >> true
            getInputs() >> [(new ProcessInput('x', String, false)): 'a', (new ProcessInput('s', Map, false)): [id: 's1']]
            getInputEnv() >> [:]
        }
        def hasher = Spy(new TaskHasher(task))
        hasher.getTaskGlobalVars() >> [:]

        expect: 'the pinned hash of a typed task without stage directives'
        hasher.compute().toString() == '6812ed701a957edc84fbe25166e015ba'
    }

    def 'should hash stage env and stdin values as strings' () {
        expect:
        stageHash(env: [FOO: 'one']) == stageHash(env: [FOO: 'one'])
        stageHash(env: [FOO: 'one']) != stageHash(env: [FOO: 'two'])
        stageHash(env: [FOO: 'one']) != stageHash(env: [BAR: 'one'])
        stageHash(env: [A: '1', B: '2']) == stageHash(env: [B: '2', A: '1'])
        and:
        stageHash(stdin: 'one') == stageHash(stdin: 'one')
        stageHash(stdin: 'one') != stageHash(stdin: 'two')
        stageKeys(stdin: 1) == stageKeys(stdin: '1')
        stageKeys(stdin: Path.of('/some/file.txt')) == ['stage_stdin', '/some/file.txt']
        and:
        stageHash(env: [X: 'one']) != stageHash(stdin: 'one')
    }

    def 'should include stage directives in the hash of typed tasks only' () {
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
        def hashFor = { boolean typed, String value ->
            def task = Mock(TaskRun) {
                getSource() >> 'hello world'
                isContainerEnabled() >> false
                getConfig() >> Mock(TaskConfig)
                getProcessor() >> processor
                hasTypedInputsOutputs() >> typed
                getInputEnv() >> [FOO: value]
            }
            def hasher = Spy(new TaskHasher(task))
            hasher.getTaskGlobalVars() >> [:]
            return hasher.compute()
        }

        expect:
        hashFor(true, 'one') == hashFor(true, 'one')
        hashFor(true, 'one') != hashFor(true, 'two')
        and: 'legacy env values are task inputs and are hashed as such'
        hashFor(false, 'one') == hashFor(false, 'two')
    }
}
