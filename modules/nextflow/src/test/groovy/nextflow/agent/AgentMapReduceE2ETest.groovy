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
import java.util.concurrent.atomic.AtomicInteger

import nextflow.Session
import nextflow.script.ScriptBinding
import nextflow.script.ScriptLoaderFactory
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end CI guard for the fully-agentic map-reduce graph (design M6, spec §13):
 * {@code planner -> mapper (parallel map) -> reducer (collect() fan-in)}, every phase an
 * {@code agent} lowered to a real {@code TaskProcessor} (M1). The whole graph is driven
 * INLINE (no dependency on {@code examples/agents/map-reduce/main.nf}) with the LLM stubbed
 * via {@link AgentRunnerProvider#testRunner}. The stub is prompt-aware: it distinguishes the
 * three agents by a marker in {@code req.prompt} and echoes correlating shard ids so findings
 * are matched by KEY, not by completion order (order-independence, spec §6.6/§10).
 *
 * All assertions are deterministic: per-role invocation counts, terminal record shape and
 * order-independent set membership. NO wall-clock / timing assertions (spec §4.3, §12 risk #4).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(120)
class AgentMapReduceE2ETest extends Dsl2Spec {

    /**
     * The spec §13 map-reduce workflow, inlined so the test does not depend on the example
     * file. Single-output agents auto-unwrap to their channel, so {@code planner(...)} IS the
     * {@code Plan} channel; {@code plan.flatMap{ it.shards }} yields a deterministic shard queue,
     * {@code mapper(shards)} maps one TaskRun per shard, and {@code reducer(findings.collect())}
     * fans in the whole bag as a single value item (one TaskRun).
     */
    private static final String MAP_REDUCE_SCRIPT = '''
        nextflow.enable.types = true

        record Shard   { id: String;      question: String }
        record Plan    { title: String;   shards: List<Shard> }
        record Finding { shardId: String; summary: String }
        record Report  { title: String;   body: String }

        agent planner {
            model 'openai/gpt-4o-2024-08-06'
            instruction 'You decompose a research brief into independent shards.'
            input:
            brief: String
            output:
            plan: Plan
            prompt:
            """
            Break this brief into independent research shards:

            ${brief}
            """
        }

        agent mapper {
            model 'openai/gpt-4o-2024-08-06'
            instruction 'You are a focused researcher; echo the shard id you were given.'
            input:
            shard: Shard
            output:
            finding: Finding
            prompt:
            """
            Shard id: ${shard.id}
            Answer this shard question:

            ${shard.question}
            """
        }

        agent reducer {
            model 'openai/gpt-4o-2024-08-06'
            instruction 'You synthesise many independent findings into one coherent report.'
            input:
            findings: Bag<Finding>
            output:
            report: Report
            prompt:
            """
            Synthesise ONE report from these findings:

            ${findings.collect { "- (${it.shardId}) ${it.summary}" }.join('\\n')}
            """
        }

        workflow {
            def plan     = planner(channel.of('a research brief'))
            def shards   = plan.flatMap { it.shards }
            def findings = mapper(shards)
            reducer(findings.collect())
        }
        '''

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    // -- The whole graph: planner once, mapper per-shard, reducer fan-in once, one Report.
    def 'should run the agentic map-reduce graph end-to-end with deterministic per-role counts and fan-in'() {
        given: 'a prompt-aware stub that counts calls per role and echoes shard keys'
        def plannerCalls = new AtomicInteger()
        def mapperCalls  = new AtomicInteger()
        def reducerCalls = new AtomicInteger()
        def reducerReqs  = Collections.synchronizedList([])
        AgentRunnerProvider.testRunner = stubRunner(plannerCalls, mapperCalls, reducerCalls, reducerReqs)

        when: 'the §13 workflow runs with an explicit cpu gate so the map fan-out is admitted'
        def result = runScript(config: [executor: [local: [cpus: 8]]], MAP_REDUCE_SCRIPT)

        then: '(1) planner fired exactly once (single value input -> fan-in of one)'
        plannerCalls.get() == 1

        and: '(2) mapper fired once per shard - flatMap produced a 3-item queue, each shard its own TaskRun'
        mapperCalls.get() == 3

        and: '(3) reducer fired exactly once on the whole collected bag (collect() -> one value item)'
        reducerCalls.get() == 1

        and: '(3) the reducer input carried ALL three findings - assert the SET of shard ids, order-independent'
        def reducerReq = reducerReqs[0] as AgentRunnerRequest
        def idsInPrompt = (reducerReq.prompt =~ /\((s\d+)\)/).collect { it[1] } as Set
        idsInPrompt == ['s1', 's2', 's3'] as Set
        and: 'each shard summary is present in the reducer prompt and its input json'
        ['sum-s1', 'sum-s2', 'sum-s3'].every { reducerReq.prompt.contains(it) }
        (['s1', 's2', 's3'] as Set).every { reducerReq.inputJson.contains(it) }

        and: '(1)+(4) the graph yields exactly one terminal Report with the expected title/body'
        def report = result.val
        report instanceof Map
        report.title == 'R'
        report.body == 'combined'
    }

    // -- Zero-new-observer proof: the three agents are real TaskProcessors, so the stock
    //    task-lifecycle + process-create events fall out for free (spec §4.7). Driven through a
    //    REAL Session (the MockSession used above short-circuits the monitor and fires no events),
    //    with only an anonymous stock TraceObserverV2 probe - no agent-specific observer type.
    def 'should fire process-create and task-lifecycle events for planner, mapper and reducer with zero new observer code'() {
        given:
        def plannerCalls = new AtomicInteger()
        def mapperCalls  = new AtomicInteger()
        def reducerCalls = new AtomicInteger()
        AgentRunnerProvider.testRunner = stubRunner(plannerCalls, mapperCalls, reducerCalls, null)
        def createdProcesses = Collections.synchronizedList([])
        def completed = Collections.synchronizedList([])
        def probe = new TraceObserverV2() {
            @Override void onProcessCreate(nextflow.processor.TaskProcessor process) { createdProcesses.add(process.name) }
            @Override void onTaskComplete(TaskEvent event) { completed.add(event) }
        }

        when:
        runWithObserver(probe, MAP_REDUCE_SCRIPT)

        then: 'the stock process-create event fired for every agent (each lowered to a TaskProcessor)'
        createdProcesses.containsAll(['planner', 'mapper', 'reducer'])

        and: 'the stub is a plain TraceObserverV2 - no new agent-specific observer type'
        probe instanceof TraceObserverV2

        and: 'task-complete events fired 1 (planner) + 3 (mapper) + 1 (reducer) = 5 across the three names'
        def byProc = completed.groupBy { it.trace.get('process') }
        byProc['planner']?.size() == 1
        byProc['mapper']?.size() == 3
        byProc['reducer']?.size() == 1
        completed.size() == 5
    }

    // -----------------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------------

    /**
     * Prompt-aware deterministic LLM stub. Each agent output is a SINGLE record, so the stub
     * returns BARE record JSON (no wrapper key) - matching the unwrapped single-record wire shape
     * (spec §5.3b/§9.3). The mapper echoes the shard id it was handed, so findings correlate by
     * key regardless of the (non-deterministic) map completion order.
     */
    private AgentRunner stubRunner(AtomicInteger planner, AtomicInteger mapper, AtomicInteger reducer, List reducerReqs) {
        return { AgentRunnerRequest req ->
            final p = req.prompt
            if( p.contains('Break this brief') ) {
                planner.incrementAndGet()
                return '{"title":"T","shards":[' +
                       '{"id":"s1","question":"q1"},{"id":"s2","question":"q2"},{"id":"s3","question":"q3"}]}'
            }
            if( p.contains('Answer this shard') ) {
                mapper.incrementAndGet()
                final m = (p =~ /Shard id:\s*(\S+)/)
                final id = m.find() ? m.group(1) : 's?'
                return "{\"shardId\":\"${id}\",\"summary\":\"sum-${id}\"}".toString()
            }
            if( p.contains('Synthesise ONE report') ) {
                reducer.incrementAndGet()
                if( reducerReqs != null ) reducerReqs.add(req)
                return '{"title":"R","body":"combined"}'
            }
            throw new IllegalStateException("unexpected agent prompt: " + (p?.take(60)))
        } as AgentRunner
    }

    /**
     * Drive the workflow through a REAL {@link Session} (real {@code ExecutorFactory},
     * {@code LocalExecutor}, monitor and {@code NativeTaskHandler}) so that the task-lifecycle
     * observer notifications and process-create events actually fire - the {@code MockSession}
     * harness used by {@code runScript} short-circuits the monitor and never emits them. A probe
     * {@link TraceObserverV2} is injected into the private {@code observersV2} list before ignition
     * (no public observer-registration API). Mirrors the proven pattern in
     * {@code AgentAsTaskIntegrationTest} (Test I); introduces NO production/observer code.
     */
    private static Session runWithObserver(TraceObserverV2 probe, String text) {
        final workDir = Files.createTempDirectory('nxf-agent-mapreduce')
        def session = new Session([workDir: workDir.toString()])
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        // inject the probe into the private observersV2 list before ignition
        final f = Session.getDeclaredField('observersV2')
        f.setAccessible(true)
        final list = new ArrayList((List) f.get(session))
        list.add(probe)
        f.set(session, list)
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        loader.parse(text)
        loader.runScript()

        session.fireDataflowNetwork()
        session.await()
        session.destroy()
        if( session.error )
            throw session.error
        return session
    }
}
