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
import java.util.concurrent.atomic.AtomicInteger

import groovy.json.JsonSlurper
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.AgentDef
import nextflow.script.BaseScript
import nextflow.script.PromptDef
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.ScriptLoaderFactory
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import spock.lang.TempDir
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockSession

import static test.ScriptHelper.runScript
import nextflow.agent.rpc.AgentRpcRegistration

/**
 * Integration tests for the tool-free agent lowered to a real {@code TaskProcessor}
 * (design M1). The LLM is stubbed via {@link AgentRunnerProvider#testRunner}; each
 * test drives a workflow through {@code runScript} (task path via the mock executor).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(30)
class AgentAsTaskIntegrationTest extends Dsl2Spec {

    @TempDir
    Path tempDir

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        Thread.interrupted() // clear any leaked interrupt flag from a fatal-tool abort (Test N)
    }

    // -- Test B: single-record output stays UNWRAPPED (bare RecordSchema.of)
    def 'should keep a single-record output unwrapped (no wrapper key)'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req; '{"answer":"ok","confidence":0.9}'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            record Answer { answer: String; confidence: Double }

            agent qa {
                model 'openai/gpt-4o'
                tools()
                input:
                q: String
                output:
                a: Answer
                prompt:
                """
                Q: ${q}
                """
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then: 'the schema is the BARE record schema - properties at the root, no wrapper key `a`'
        captured.outputSchema.type == 'object'
        captured.outputSchema.properties.containsKey('answer')
        captured.outputSchema.properties.containsKey('confidence')
        !captured.outputSchema.properties.containsKey('a')
        and: 'the bare record schema carries required + additionalProperties:false (byte shape unchanged)'
        (captured.outputSchema.required as Set) == ['answer', 'confidence'] as Set
        captured.outputSchema.additionalProperties == false
        and: 'inputJson is the bare single value (not a {q:...} wrapper)'
        captured.inputJson == '"hello"'
        and: 'the channel yields the parsed record'
        def out = result.val
        out instanceof Map
        out.answer == 'ok'
        out.confidence == 0.9d
    }

    // -- Test C: free-text passthrough (single scalar output)
    def 'should pass through free text for a single scalar output'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req; 'hello world'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                tools()
                input:
                q: String
                output:
                answer: String
                prompt:
                """
                Q: ${q}
                """
            }

            workflow {
                qa(channel.of('hi'))
            }
            ''')

        then:
        result.val == 'hello world'
        captured.outputSchema == null
    }

    // -- Test D: multiple inputs combine natively
    def 'should combine multiple inputs and render both in the prompt + inputJson map'() {
        given:
        def prompts = Collections.synchronizedList([])
        def jsons = Collections.synchronizedList([])
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            prompts.add(req.prompt); jsons.add(req.inputJson); 'ok'
        } as AgentRunner

        when:
        runScript('''
            nextflow.enable.types = true

            agent combiner {
                model 'openai/gpt-4o'
                tools()
                input:
                a: String
                b: String
                output:
                r: String
                prompt:
                """
                A=${a} B=${b}
                """
            }

            workflow {
                combiner(channel.of('x1','x2'), channel.of('y1','y2'))
            }
            ''')

        then: 'one invocation per position (2), each rendering both inputs'
        prompts.size() == 2
        prompts.every { it.contains('A=') && it.contains('B=') }
        and: 'for N>1 inputs the inputJson is a {a:..,b:..} map'
        jsons.every { it.contains('"a"') && it.contains('"b"') }
        (prompts.join(' ')).contains('x1')
        (prompts.join(' ')).contains('y1')
    }

    // -- Test E1: queue input runs per-item
    def 'should run once per item for a queue input'() {
        given:
        def count = new AtomicInteger()
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            count.incrementAndGet(); "r-${req.inputJson}".toString()
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                tools()
                input:
                q: String
                output:
                answer: String
                prompt: "Q: ${q}"
            }

            workflow {
                qa(channel.of('x','y','z'))
            }
            ''')

        then:
        count.get() == 3
        and:
        def ch = result
        def vals = [ch.val, ch.val, ch.val]
        vals.size() == 3
    }

    // -- Test E2 [FAN-IN]: value/singleton input (collect()) runs exactly once
    def 'should run exactly once for a value/singleton input (fan-in)'() {
        given:
        def count = new AtomicInteger()
        Object capturedInput = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            count.incrementAndGet(); capturedInput = req.inputJson; 'summary'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent reducer {
                model 'openai/gpt-4o'
                tools()
                input:
                items: String
                output:
                report: String
                prompt: "Reduce: ${items}"
            }

            workflow {
                reducer(channel.of('a','b','c').collect())
            }
            ''')

        then: 'the agent fired exactly once with the whole collected bag'
        count.get() == 1
        capturedInput.contains('a') && capturedInput.contains('b') && capturedInput.contains('c')
        and:
        result.val == 'summary'
    }

    // -- Test I: parity smoke - the agent runs as a genuine TaskProcessor
    def 'should run as a real TaskProcessor firing process-create and task-lifecycle events with a per-task work dir (zero new observer code)'() {
        given:
        def ran = new java.util.concurrent.atomic.AtomicBoolean(false)
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> ran.set(true); 'done' } as AgentRunner
        def createdProcesses = Collections.synchronizedList([])
        def submitted = Collections.synchronizedList([])
        def completed = Collections.synchronizedList([])
        def probe = new TraceObserverV2() {
            @Override void onProcessCreate(TaskProcessor process) { createdProcesses.add(process.name) }
            @Override void onTaskSubmit(TaskEvent event) { submitted.add(event) }
            @Override void onTaskComplete(TaskEvent event) { completed.add(event) }
        }

        when:
        def session = runWithObserver(probe, '''
            nextflow.enable.types = true

            agent probe_agent {
                model 'openai/gpt-4o'
                tools()
                input:
                q: String
                output:
                answer: String
                prompt: "Q: ${q}"
            }

            workflow {
                probe_agent(channel.of('hi'))
            }
            ''')

        then: 'the agent lowered to a TaskProcessor that emitted the standard process-create event'
        createdProcesses.contains('probe_agent')
        and: 'the task body actually executed (LLM call ran in-JVM)'
        ran.get()
        and: 'task-lifecycle events fired for the agent-named task (falls out for free from the monitor)'
        def submit = submitted.find { it.trace.get('process') == 'probe_agent' }
        submit != null
        def complete = completed.find { it.trace.get('process') == 'probe_agent' }
        complete != null
        and: 'the handler carries a per-task work dir under the session work dir'
        def workDir = complete.handler.task.workDir
        workDir != null
        Files.exists(workDir)
        workDir.startsWith(session.workDir)
    }

    // -- Test J: parallel map with explicit executor.local.cpus (count/independence only)
    def 'should map over a queue producing one output per item (parallelism config, no timing assertion)'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            "out-${req.inputJson}".toString()
        } as AgentRunner

        when:
        def result = runScript(
            config: [executor: [local: [cpus: 4]]],
            '''
            nextflow.enable.types = true

            agent mapper {
                model 'openai/gpt-4o'
                tools()
                input:
                n: Integer
                output:
                r: String
                prompt: "N: ${n}"
            }

            workflow {
                mapper(channel.of(1,2,3,4,5,6,7,8))
            }
            ''')

        then: 'all 8 tasks complete with 8 distinct outputs'
        def vals = (1..8).collect { result.val }
        vals.size() == 8
        (vals as Set).size() == 8
    }

    // -- Test K [M-Skills]: a skills-ONLY agent (declares `skills`, NO `tools`) lowers to the
    //    TASK path - proven by the process-create event firing for the agent-named task and the
    //    map fan-out invoking the runner once per queue item. On CURRENT source a skills agent is
    //    grouped with tools onto the legacy operator path, so NO TaskProcessor is created and
    //    `createdProcesses` never contains the agent name (RED).
    def 'should run a skills-only agent on the task path (process-create + map fan-out)'() {
        given: 'a local skill fixture beside the script (skills/<name>/SKILL.md with YAML frontmatter)'
        def skillDir = Files.createDirectories(tempDir.resolve('skills').resolve('classifier'))
        Files.writeString(skillDir.resolve('SKILL.md'),
            "---\nname: classifier\ndescription: labels an item\n---\nClassify the given item into a label.")
        def script = tempDir.resolve('main.nf')
        Files.writeString(script, '''
            nextflow.enable.types = true

            agent classify {
                model 'openai/gpt-4o'
                skills 'classifier'
                input:
                item: String
                output:
                label: String
                prompt: "Classify: ${item}"
            }

            workflow {
                classify(channel.of('a','b'))
            }
            ''')

        and: 'a runner stub that records how many times it fires and what skills each request carries'
        def calls = new AtomicInteger()
        def capturedSkills = Collections.synchronizedList([])
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            calls.incrementAndGet()
            capturedSkills.add(req.skills)
            'ok'
        } as AgentRunner

        and: 'a probe observer to capture the process-create event (task-path-only signal)'
        def createdProcesses = Collections.synchronizedList([])
        def probe = new TraceObserverV2() {
            @Override void onProcessCreate(TaskProcessor process) { createdProcesses.add(process.name) }
        }

        when:
        runWithObserver(probe, script)

        then: 'the skills-only agent lowered to a real TaskProcessor (NOT the legacy operator)'
        createdProcesses.contains('classify')
        and: 'the map fanned out: one runner invocation per queue item (2 inputs -> 2 calls)'
        calls.get() == 2
        and: 'every request carried the resolved skill descriptor'
        capturedSkills.size() == 2
        capturedSkills.every { it != null && (it*.name as Set) == ['classifier'] as Set }
    }

    // -- Test L [M-Tools (a)]: a TOOLS agent (declares `tools 'nf:module_run:greet'`, NO skills) now RUNS ON THE
    //    TASK path. Proven by the process-create event firing for the agent-named task, the map
    //    fanning out (2 inputs -> 2 runner invocations), the tool call round-tripping through the
    //    REAL bridge (req.toolSpecs/req.dispatch populated -> drives the real `greet` process), AND
    //    session.await() COMPLETING (reaching the `then` block proves the bridge was closed and the
    //    run did not hang). On CURRENT source a tools agent is routed to the legacy GPars operator
    //    (runLegacy), so NO TaskProcessor is created and `createdProcesses` never contains the agent
    //    name (RED).
    def 'should run a tools agent on the task path (process-create + map fan-out + tool round-trip + no hang)'() {
        given: 'a runner stub that drives the REAL bridge and records the fan-out'
        def calls = new AtomicInteger()
        def dispatchPresent = Collections.synchronizedList([])
        def toolNames = Collections.synchronizedList([])
        def dispatchResults = Collections.synchronizedList([])
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            calls.incrementAndGet()
            dispatchPresent.add(req.dispatch != null)
            toolNames.addAll(req.toolSpecs*.name)
            final name = new JsonSlurper().parseText(req.inputJson).toString()
            // Each agent request invokes the same tool. Request-scoped channels keep the
            // reversed/concurrent completions correlated with the correct caller.
            final result = req.dispatch.call('greet', """{"name":"${name}"}""")
            dispatchResults.add([name: name, result: result])
            'answer'
        } as AgentRunner

        and: 'a probe observer to capture the process-create event (task-path-only signal)'
        def createdProcesses = Collections.synchronizedList([])
        def activeTools = new AtomicInteger()
        def peakTools = new AtomicInteger()
        def probe = new TraceObserverV2() {
            @Override void onProcessCreate(TaskProcessor process) { createdProcesses.add(process.name) }
            @Override void onTaskStart(TaskEvent event) {
                if( event.handler.task.processor.name == 'greet' ) {
                    final count = activeTools.incrementAndGet()
                    peakTools.updateAndGet { int peak -> Math.max(peak, count) }
                }
            }
            @Override void onTaskComplete(TaskEvent event) {
                if( event.handler.task.processor.name == 'greet' )
                    activeTools.decrementAndGet()
            }
        }

        when:
        runWithObserver(probe, '''
            nextflow.enable.types = true

            process greet {
                input:
                name: String

                output:
                greeting: String

                exec:
                Thread.sleep(name == 'a' ? 300 : 100)
                greeting = "Hello ${name}!"
            }

            agent assistant {
                model 'openai/gpt-4o'
                tools 'nf:module_run:greet'
                input:
                request: String
                output:
                answer: String
                prompt: "Handle: ${request}"
            }

            workflow {
                assistant(channel.of('a','b'))
            }
            ''')

        then: 'the tools agent lowered to a real TaskProcessor (NOT the legacy operator)'
        createdProcesses.contains('assistant')
        and: 'the map fanned out: one runner invocation per queue item (2 inputs -> 2 calls)'
        calls.get() == 2
        and: 'each request carried the live tool descriptors + dispatch callback'
        dispatchPresent.every { it }
        (toolNames as Set) == ['greet'] as Set
        and: 'same-tool calls executed concurrently and each reply stayed correlated'
        peakTools.get() == 2
        dispatchResults.collectEntries {
            [(it.name): new JsonSlurper().parseText(it.result).greeting]
        } == [a: 'Hello a!', b: 'Hello b!']
        // reaching this point proves session.await() COMPLETED (the bridge was poisoned/closed);
        // a leaked, un-poisoned tool queue would hang await() and trip the class @Timeout instead.
    }

    // -- Test M [M-Tools (b)]: agents are cacheable, with or without tools. Uses the white-box
    //    buildAgentTask path (the same seam AgentResumeIntegrationTest drives) to assert the
    //    isCacheable() gate directly.
    def 'an agent processor is cacheable with tools, without tools, and under a launch spec'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner
        newSession()

        when: 'an fs:-tool agent (the family needs no in-scope process) is lowered'
        def toolsProc = newAgent([model: 'openai/gpt-4o', tools: ['fs:*']]).buildAgentTask(['hello'])

        then:
        toolsProc.getConfig().isCacheable() == true

        when: 'a tool-free agent is lowered'
        def freeProc = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'the tool-free agent keeps the default cacheable behaviour (unchanged)'
        freeProc.getConfig().isCacheable() == true

        when: 'a tool-free external runner with a canonical launch spec is lowered'
        AgentRunnerProvider.testRunner = new AgentRunner() {
            @Override
            String getName() { 'external' }

            @Override
            AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(
                    containerProxyCommand: ['/opt/proxy'],
                    containerHarnessCommand: ['/opt/harness'])
            }

            @Override
            String run(AgentRunnerRequest request) { 'x' }
        }
        newSession(containerized())
        def externalProc = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'using a launch spec does not disable resume for a tool-free agent'
        externalProc.getConfig().isCacheable() == true
    }

    // -- Test M3 [canonical launch path]: the broker lives in the runner plugin, so core asks the
    //    RESOLVED runner for the registration (AgentRunner.register) and splices it into the launch
    //    command at the `--` separator. Drives the task body directly with a stub runner: the whole
    //    core-side contract is the registration call plus the flag names, and neither is covered by
    //    the plugin-side broker test.
    def 'should ask the resolved runner to register and splice the result into the launch command'() {
        given:
        def registered = []
        AgentRunnerProvider.testRunner = canonicalRunner(registered)
        newSession(containerized())

        when: 'a canonical agent task body runs'
        def script = runTaskBody(newAgent([model: 'openai/gpt-4o']), new TaskConfig([container: 'agent-image:test']))

        then: 'a canonical task is ALWAYS containerized, so the endpoint is ALWAYS dialled remotely'
        registered == [[prompt: 'Q', remote: true]]

        and: 'the connection flags are spliced before the `--`, leaving the harness command last'
        // the pinned certificate digest must be there: it is what stops a silent regression to an
        // unpinned, cleartext connection, and unlike the token it is public so the script is fine
        script == "exec '/opt/agent-rpc' '--log' 'debug' " +
            "'--endpoint' 'host.docker.internal:9999' '--invocation' 'inv-1' '--fingerprint' 'abc123' '--token' 'tok-1' " +
            "'--' 'node' '/opt/runner.mjs'"
    }

    // -- the build-time containerization guard cannot see a LAZY `agent.container` (a closure or a
    //    GString resolving per task), so the body re-checks the resolved image. Without this the
    //    in-image proxy path would be exec'd on the host as `No such file or directory`.
    def 'a canonical task whose container resolves to nothing fails loudly, never exec-ing host paths'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner([])
        newSession(containerized())

        when: 'the per-task container resolves to null even though the build-time value was truthy'
        runTaskBody(newAgent([model: 'openai/gpt-4o']), new TaskConfig([:]))

        then:
        def err = thrown(nextflow.exception.ScriptRuntimeException)
        err.message.contains('agent.container')
    }

    // -- an agent whose configuration would NOT run the task in a container is rejected BEFORE the
    //    run starts: the launch command is made of in-image paths only.
    def 'a canonical agent is rejected when the resolved configuration would not containerize it'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner([])

        when: 'no `agent.container` is declared'
        newSession([docker: [enabled: true]])
        newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then:
        def missing = thrown(nextflow.exception.ScriptRuntimeException)
        missing.message.contains('must declare a container')

        when: 'an image is declared but no container engine is enabled'
        newSession([agent: [container: 'agent-image:test']])
        newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then:
        def disabled = thrown(nextflow.exception.ScriptRuntimeException)
        disabled.message.contains('would not run it in a container')
    }

    // -- the agent task dials the driver's broker back at `host.docker.internal` when no explicit
    //    address is configured; on Linux Docker that name resolves ONLY when the container is run
    //    with `--add-host=...:host-gateway`, so the driver adds it to the agent task itself instead
    //    of leaving the most common local setup to fail with a connection timeout.
    def 'a locally containerized agent task is given the docker host-gateway run option'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner([])

        when: 'docker runs the container and no `agent.rpc.remoteHost` is set'
        newSession(containerized())
        def proc = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then:
        proc.getConfig().get('containerOptions') == '--add-host=host.docker.internal:host-gateway'

        when: 'an explicit driver address is configured instead'
        newSession(containerized([agent: [rpc: [remoteHost: '127.0.0.1']]]))
        proc = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'nothing is added - the alias is not what the task uses'
        proc.getConfig().get('containerOptions') == null

        when: 'the agent declares container options of its own'
        newSession(containerized([agent: [containerOptions: '--cpus 2']]))
        proc = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'they are preserved, never replaced'
        proc.getConfig().get('containerOptions') == '--cpus 2 --add-host=host.docker.internal:host-gateway'
    }

    // -- Test M4: a runner that publishes a launch spec but no register() override fails loudly
    //    (the AgentRunner default), instead of silently launching a proxy with no endpoint.
    def 'should fail with a clear error when a launch-spec runner does not implement register'() {
        given:
        AgentRunnerProvider.testRunner = new AgentRunner() {
            @Override
            String getName() { 'no-broker' }

            @Override
            AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(containerProxyCommand: ['/opt/proxy'], containerHarnessCommand: ['/opt/harness'])
            }

            @Override
            String run(AgentRunnerRequest request) { 'x' }
        }
        newSession(containerized())

        when:
        runTaskBody(newAgent([model: 'openai/gpt-4o']), new TaskConfig([container: 'agent-image:test']))

        then:
        def err = thrown(UnsupportedOperationException)
        err.message == 'Agent runner `no-broker` provides a launch spec but no RPC broker'
    }

    // -- Test M2: process selectors do not configure agents.
    // -- an agent is admitted without a cpu/capacity throttle, so absent a cap it fans out as wide
    //    as its input channel. A default bounds concurrent LLM calls; the user's value still wins.
    def 'an agent gets a default maxForks that an explicit value overrides'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner

        when: 'nothing is declared'
        newSession()
        def defaulted = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then:
        defaulted.getMaxForks() == AgentConfig.DEFAULT_MAX_FORKS

        when: 'the agent scope raises it'
        newSession([agent: [maxForks: 42]])
        def raised = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'the user value wins over the default'
        raised.getMaxForks() == 42
    }

    def 'an agent ignores a process selector that sets cache false'() {
        given: 'a process selector that would have matched the agent under the old shared scope'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner
        newSession([process: ['withName:qa': [cache: false]]])

        when: 'a tools agent named qa is lowered from the independent agent scope'
        def toolsProc = newAgent([model: 'openai/gpt-4o', tools: ['fs:*']]).buildAgentTask(['hello'])

        then: 'the process selector is irrelevant and the agent keeps its default cacheability'
        toolsProc.getConfig().isCacheable() == true
    }

    // -- the counterpart: the same setting in the `agent` scope IS honoured.
    def 'an agent honours an agent selector that sets cache false'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner
        newSession([agent: ['withName:qa': [cache: false]]])

        when:
        def toolsProc = newAgent([model: 'openai/gpt-4o', tools: ['fs:*']]).buildAgentTask(['hello'])

        then:
        toolsProc.getConfig().isCacheable() == false
    }

    // -- the RESOLVED Executor is the oracle for agent placement, not the config map: the
    //    TaskProcessor holds the Executor instance and never re-derives it from the config, so
    //    an `executor` write that lands AFTER createTaskProcessorResolved silently does nothing.
    //    Asserting only config.executor would keep passing if the ordering regressed, while every
    //    in-JVM agent quietly ran on the compute executor -- and an agent body blocks on its own
    //    tool sub-tasks, so that deadlocks concurrent tool agents.
    def 'an in-JVM agent is resolved onto the agent executor, not merely configured for it'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner
        newSession()

        when:
        def proc = newAgent([model: 'openai/gpt-4o', tools: ['fs:*']]).buildAgentTask(['hello'])

        then: 'the Executor the processor actually holds'
        proc.getExecutor().getName() == 'agent'
        and: 'and the config it was resolved from'
        proc.getConfig().get('executor') == 'agent'

        when: 'a launch-spec runner is used instead -- it is NOT pinned to the agent executor'
        // a FRESH session is required: ExecutorFactory caches the resolved Executor by class and
        // MockExecutorFactory maps every executor name to MockExecutor, so reusing the session
        // above would hand back the instance already named `agent` and the check would be vacuous
        newSession(containerized())
        AgentRunnerProvider.testRunner = new AgentRunner() {
            @Override
            String getName() { 'external' }

            @Override
            AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(
                    containerProxyCommand: ['/opt/proxy'],
                    containerHarnessCommand: ['/opt/harness'])
            }

            @Override
            String run(AgentRunnerRequest request) { 'x' }
        }
        def offloadable = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])

        then: 'it resolves to the standard agent default, proving the assertion above discriminates'
        offloadable.getExecutor().getName() == AgentConfig.DEFAULT_EXECUTOR
    }

    // -- Test N [M-Tools (c)]: a hard-failing tool aborts the run cleanly (the AgentToolFatalError
    //    path) rather than hanging. Deterministically reproduces the session-abort interrupt by
    //    swapping the wired tool's captured output channel with a stub whose getVal() throws the
    //    GPars-wrapped InterruptedException (== a session abort on task failure), exactly as
    //    ModuleToolBridgeTaskFailureTest does. Here the AgentToolFatalError is left UNCAUGHT so it
    //    propagates out of the runner and must abort the run. On CURRENT source the agent runs on
    //    the legacy operator, so the process-create event never fires for the agent (RED); the run
    //    must still abort (not hang), surfacing an AgentToolFatalError.
    def 'a hard-failing request-scoped tool aborts the run cleanly on the task path'() {
        given: 'a runner stub that swaps the tool output for one that throws the abort interrupt, then calls the tool (uncaught)'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            // do NOT catch: the failed process aborts the session, interrupts the
            // correlated reply pull and must surface as AgentToolFatalError.
            req.dispatch.call('greet', '{"name":"Ada"}')
        } as AgentRunner

        and: 'a probe observer to capture the process-create event (task-path-only signal)'
        def createdProcesses = Collections.synchronizedList([])
        def probe = new TraceObserverV2() {
            @Override void onProcessCreate(TaskProcessor process) { createdProcesses.add(process.name) }
        }

        when:
        Throwable err = null
        try {
            runWithObserver(probe, '''
                nextflow.enable.types = true

                process greet {
                    input:
                    name: String

                    output:
                    greeting: String

                    exec:
                    throw new IllegalStateException('tool failed')
                }

                agent assistant {
                    model 'openai/gpt-4o'
                    tools 'nf:module_run:greet'
                    input:
                    request: String
                    output:
                    answer: String
                    prompt: "Handle: ${request}"
                }

                workflow {
                    assistant(channel.of('hi'))
                }
                ''')
        }
        catch( Throwable t ) {
            err = t
        }

        then: 'the tools agent lowered to a real TaskProcessor (NOT the legacy operator)'
        createdProcesses.contains('assistant')
        createdProcesses.contains('greet')
        and: 'the fatal tool error aborted the run (it did not hang and did not silently succeed)'
        err != null
    }

    // -- Test P: endpoint/credential plumbing (design D4). The ladder is resolved ONCE, in core,
    //    and handed to whichever runner is selected -- no runner reads the environment.
    def 'should hand the resolved endpoint and credential to the runner without leaking the key'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'ok' } as AgentRunner
        and: 'an empty environment, so only the config can supply a value'
        SysEnv.push([:])

        when:
        def result = runScript(config: [agent: [apiKey: 'sk-canary-51ad', baseUrl: 'http://localhost:8000/v1']], '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                tools()
                input:
                q: String
                output:
                answer: String
                prompt: "Q: ${q}"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then: 'both reach the runner through the request'
        result.val == 'ok'
        captured.apiKey == 'sk-canary-51ad'
        captured.baseUrl == 'http://localhost:8000/v1'

        and: 'but the portable payload that crosses the plaintext RPC link carries only the endpoint'
        AgentProtocolSpec.fromRequest(captured).baseUrl == 'http://localhost:8000/v1'
        !AgentProtocolSpec.fromRequest(captured).containsKey('apiKey')

        and: 'and an interpolated request cannot leak the credential into the log'
        !"$captured".toString().contains('sk-canary-51ad')

        cleanup:
        SysEnv.pop()
    }

    // -----------------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------------

    /**
     * The minimal configuration that CONTAINERIZES a canonical agent task -- an enabled container
     * engine plus an `agent.container` image -- which every launch-spec runner now requires, since
     * its launch command is built from paths that exist only inside the runner image.
     */
    private static Map containerized(Map config = [:]) {
        final result = new LinkedHashMap(config)
        result.docker = [enabled: true]
        final agentScope = new LinkedHashMap((Map) (config.agent ?: [:]))
        agentScope.container = 'agent-image:test'
        result.agent = agentScope
        return result
    }

    /**
     * A launch-spec runner publishing only in-image paths, recording each registration into
     * {@code registered} so a test can assert what core asked the broker for.
     */
    private static AgentRunner canonicalRunner(List registered) {
        return new AgentRunner() {
            @Override
            String getName() { 'external' }

            @Override
            AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(
                    containerProxyCommand: ['/opt/agent-rpc', '--log', 'debug'],
                    containerHarnessCommand: ['node', '/opt/runner.mjs'])
            }

            @Override
            AgentRpcRegistration register(AgentRunnerRequest request, boolean remote) {
                registered << [prompt: request.prompt, remote: remote]
                // a fingerprint (or the explicit --insecure opt-out) is mandatory: transportArgs()
                // rejects a registration carrying neither rather than dialling unpinned
                return new AgentRpcRegistration('inv-1', 'tok-1', 'host.docker.internal:9999', 'abc123')
            }

            @Override
            String run(AgentRunnerRequest request) { throw new UnsupportedOperationException('canonical task path') }
        }
    }

    /** Spin a MockSession (MockExecutorFactory) and make it the global session (white-box build path). */
    private Session newSession(Map config = null) {
        def session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()
        Global.session = session
        return session
    }

    /**
     * Lower the agent and invoke the synthesized task body against a stand-in task context, so the
     * canonical launch command can be asserted without igniting the dataflow network. The body is
     * DELEGATE_ONLY over the task context: a plain map supplying the declared input and `task` is
     * all it reads.
     */
    private static String runTaskBody(AgentDef agent, TaskConfig taskConfig) {
        final body = (Closure) agent.buildAgentTask(['hello']).getTaskBody().closure.clone()
        body.setDelegate([q: 'hello', task: taskConfig])
        body.setResolveStrategy(Closure.DELEGATE_ONLY)
        return body.call()
    }

    /** Construct an {@link AgentDef} directly (single String in/out, trivial prompt) for white-box builds. */
    private AgentDef newAgent(Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    /**
     * Path-based variant of {@link #runWithObserver(TraceObserverV2, String)}: initializes the real
     * {@link Session} from a {@link ScriptFile} so the script's base/module dir is set and a local
     * {@code skills/} directory beside the script resolves. Injects the probe observer into the
     * private {@code observersV2} list before ignition (no public registration API), then runs.
     */
    private static Session runWithObserver(TraceObserverV2 probe, Path script) {
        final workDir = Files.createTempDirectory('nxf-agent-skills-test')
        def session = new Session([workDir: workDir.toString()])
        session.setBinding(new ScriptBinding())
        session.init(new ScriptFile(script), null, null, null)
        final f = Session.getDeclaredField('observersV2')
        f.setAccessible(true)
        final list = new ArrayList((List) f.get(session))
        list.add(probe)
        f.set(session, list)
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        loader.parse(script)
        loader.runScript()

        session.fireDataflowNetwork()
        session.await()
        session.destroy()
        if( session.error )
            throw session.error
        return session
    }

    /** Concatenate the message of a throwable and its cause chain. */
    private static String allMessages(Throwable t) {
        final sb = new StringBuilder()
        while( t != null ) {
            if( t.message ) sb.append(t.message).append(' | ')
            t = t.cause
        }
        return sb.toString()
    }

    /**
     * Drive the workflow through a *real* {@link Session} (real {@code ExecutorFactory},
     * {@code LocalExecutor}, {@code TaskPollingMonitor} and {@code NativeTaskHandler}) so that
     * the task-lifecycle observer notifications ({@code notifyTaskSubmit/Start/Complete}) and a
     * real per-task work dir actually fall out of the run — the MockSession harness used by the
     * other tests short-circuits the monitor and never fires those events. A probe
     * {@link TraceObserverV2} is injected into the private {@code observersV2} list before
     * ignition (there is no public observer-registration API). Uses an isolated temp work dir
     * and returns the live {@link Session} so callers can assert on {@code session.workDir}.
     */
    private static Session runWithObserver(TraceObserverV2 probe, String text) {
        final workDir = Files.createTempDirectory('nxf-agent-test')
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
