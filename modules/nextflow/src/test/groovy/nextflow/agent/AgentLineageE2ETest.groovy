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
package nextflow.agent

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import nextflow.Session
import nextflow.lineage.DefaultLinStore
import nextflow.lineage.LinObserver
import nextflow.lineage.LinPropertyValidator
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.model.v1beta1.AgentRun
import nextflow.lineage.model.v1beta1.TaskOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.ScriptLoaderFactory
import spock.lang.TempDir
import spock.lang.Timeout
import test.Dsl2Spec

/**
 * End-to-end test of lineage capture for the {@code agent} primitive: runs a real workflow
 * through a real {@link Session} (real executor, monitor and task handlers, so the task
 * lifecycle notifications actually fire) with a {@link LinObserver} writing to a local
 * lineage store, then asserts on what landed on disk.
 *
 * <p>The LLM is stubbed via {@link AgentRunnerProvider#testRunner} — this test is in package
 * {@code nextflow.agent} so that package-scoped assignment is a plain field write.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(120)
class AgentLineageE2ETest extends Dsl2Spec {

    @TempDir
    Path tempDir

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        AgentCallInfo.clear()
    }

    def 'should record an agent run as an AgentRun lineage record'() {
        given: 'a stubbed LLM that also reports a concrete resolved model'
        AgentRunnerProvider.testRunner = new AgentRunner() {
            @Override String getName() { 'stub-runner' }
            @Override String run(AgentRunnerRequest req) {
                AgentCallInfo.setResolvedModel('gpt-4o-2024-11-20')
                return 'the answer is 42'
            }
        }

        and: 'a workflow with one agent and one ordinary process'
        final script = tempDir.resolve('main.nf')
        script.text = '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                instruction 'You are terse.'
                maxIterations 7
                input:
                q: String
                output:
                answer: String
                prompt:
                """
                Question: ${q}
                """
            }

            process plain {
                input:
                x: String

                output:
                stdout()

                script:
                """
                echo ${x}
                """
            }

            workflow {
                qa(channel.of('what is the meaning of life'))
                plain(channel.of('hello'))
            }
            '''

        when:
        final store = runWithLineage(script)

        then: 'the agent task is stored as an AgentRun, not a TaskRun'
        final agentRun = loadOnly(store, AgentRun)
        agentRun.name.startsWith('qa')
        agentRun.sessionId
        agentRun.codeChecksum

        and: 'it records the resolved agent identity'
        agentRun.runner == 'stub-runner'
        agentRun.model == 'openai/gpt-4o'
        agentRun.instruction == 'You are terse.'
        agentRun.maxIterations == 7
        agentRun.promptTemplate.contains('Question:')
        agentRun.tools == null
        agentRun.skills == null
        agentRun.workflowRun.startsWith('lid://')

        and: 'the concrete model reported by the provider is captured'
        agentRun.resolvedModel == 'gpt-4o-2024-11-20'

        and: 'the typed input is recorded as a value parameter, not as the literal class name'
        agentRun.input.size() == 1
        agentRun.input[0].name == 'q'
        agentRun.input[0].type == 'val'
        agentRun.input[0].value == 'what is the meaning of life'

        and: 'the agent output is recorded as a TaskOutput, exactly like a process'
        final outputs = loadAll(store, TaskOutput)
        final agentOutput = outputs.find { it.taskRun == "lid://${agentHash(store)}".toString() }
        agentOutput
        agentOutput.output[0].name == 'answer'
        agentOutput.output[0].value == 'the answer is 42'

        and: 'the ordinary process in the same run is still stored as a plain TaskRun'
        final taskRuns = loadAll(store, TaskRun)
        taskRuns.size() == 1
        taskRuns[0].name.startsWith('plain')
        taskRuns[0].script.contains('echo')

        and: 'the record is queryable by its agent-specific fields, i.e. `lineage find type=AgentRun model=...`'
        new LinPropertyValidator().validateQueryParams(['type', 'model', 'runner'] as Set)
        store.search([type: ['AgentRun'], model: ['openai/gpt-4o']]).toList() == [agentHash(store)]

        and: 'the stored agent record carries no script field at all -- on the RPC runner path a\n            rendered agent command embeds a per-invocation capability token'
        final json = rawJson(store, agentHash(store))
        json.kind == 'AgentRun'
        json.version == 'lineage/v1beta1'
        !json.spec.containsKey('script')
    }

    def 'should record the tools and skills an agent was allowed to use'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            // report a snapshot: asserting it is NOT recorded below then proves the
            // non-cacheable suppression, rather than passing because nothing set it
            AgentCallInfo.setResolvedModel('gpt-4o-2024-11-20')
            return 'done'
        } as AgentRunner

        and:
        final script = tempDir.resolve('main.nf')
        script.text = '''
            nextflow.enable.types = true

            process upper {
                input:
                s: String

                output:
                shouted: String

                exec:
                shouted = s.toUpperCase()
            }

            agent shouty {
                model 'openai/gpt-4o'
                tools 'nf:module_run:upper'
                input:
                word: String
                output:
                result: String
                prompt:
                """
                Shout: ${word}
                """
            }

            workflow {
                shouty(channel.of('hi'))
            }
            '''

        when:
        final store = runWithLineage(script)

        then:
        final agentRun = loadOnly(store, AgentRun)
        agentRun.tools == ['upper']
        agentRun.skills == null
        // a module-tool agent IS cacheable (its tools are folded into the cache key by
        // AgentDef.toolsFingerprint), so the snapshot the runner reported is persisted for
        // drift detection just like a tool-free agent's
        agentRun.model == 'openai/gpt-4o'
        agentRun.resolvedModel == 'gpt-4o-2024-11-20'
    }

    def 'should not let a process selector override the agent marker'() {
        // createTaskProcessor applies the process config scope onto the SAME config object, so
        // without the re-assert in buildAgentTask a selector would silently downgrade the record
        given: 'a config that tries to replace the agent identity through a withName: selector'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner

        and:
        final script = tempDir.resolve('main.nf')
        script.text = '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                answer: String
                prompt: "Q: ${q}"
            }

            workflow {
                qa(channel.of('hello'))
            }
            '''

        when:
        final store = runWithLineage(script, [process: ['withName:qa': [agentInfo: 'forged']]])

        then: 'the run is still recorded as a genuine AgentRun, not as a plain TaskRun'
        final agentRun = loadOnly(store, AgentRun)
        agentRun.model == 'openai/gpt-4o'
        and: 'and the forged value never reaches the store'
        loadAll(store, TaskRun).isEmpty()
    }

    // ---- helpers -------------------------------------------------------------------------

    /**
     * Run the script through a real {@link Session} with a {@link LinObserver} writing to a
     * local lineage store under the temp dir, and return the opened store.
     *
     * <p>The observer is constructed directly and injected into the private {@code observersV2}
     * list rather than discovered through the plugin system: this keeps the test independent of
     * plugin-manager lifecycle state leaking in from other tests in the same JVM.
     */
    private DefaultLinStore runWithLineage(Path script, Map extraConfig = [:]) {
        final workDir = Files.createTempDirectory('nxf-agent-lineage-test')
        final storeDir = tempDir.resolve('lineage')
        final config = [
            workDir: workDir.toString(),
            lineage: [enabled: true, store: [location: storeDir.toString()]] ] + extraConfig

        final session = new Session(config)
        session.setBinding(new ScriptBinding())
        // MUST init from a real ScriptFile: LinObserver.collectScriptDataPaths dereferences
        // workflowMetadata.scriptFile unguarded, which is null when initialized from null
        session.init(new ScriptFile(script), null, null, null)

        final store = new DefaultLinStore()
        store.open(LineageConfig.create(session))
        injectObserver(session, new LinObserver(session, store))

        session.start()
        final loader = ScriptLoaderFactory.create(session)
        loader.parse(script)
        loader.runScript()
        session.fireDataflowNetwork()
        session.await()
        session.destroy()
        if( session.error )
            throw session.error
        return store
    }

    private static void injectObserver(Session session, Object probe) {
        final f = Session.getDeclaredField('observersV2')
        f.setAccessible(true)
        final list = new ArrayList((List) f.get(session))
        list.add(probe)
        f.set(session, list)
    }

    /** The store keys are the record directories holding a `.data.json` file. */
    private List<String> recordKeys(DefaultLinStore store) {
        return Files.list(tempDir.resolve('lineage'))
            .filter { Files.exists(it.resolve('.data.json')) }
            .map { it.fileName.toString() }
            .sorted()
            .toList()
    }

    private <T> List<T> loadAll(DefaultLinStore store, Class<T> type) {
        return recordKeys(store).collect { store.load(it) }.findAll { type.isInstance(it) } as List<T>
    }

    private <T> T loadOnly(DefaultLinStore store, Class<T> type) {
        final found = loadAll(store, type)
        assert found.size() == 1, "Expected exactly one ${type.simpleName} record, found ${found.size()}"
        return found[0]
    }

    private String agentHash(DefaultLinStore store) {
        return recordKeys(store).find { store.load(it) instanceof AgentRun }
    }

    private Map rawJson(DefaultLinStore store, String key) {
        return new JsonSlurper().parse(tempDir.resolve("lineage/${key}/.data.json")) as Map
    }
}
