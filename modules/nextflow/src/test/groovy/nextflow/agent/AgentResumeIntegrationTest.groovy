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

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskStartParams
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.AgentDef
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.PromptDef
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceRecord
import nextflow.trace.event.TaskEvent
import nextflow.util.CacheHelper
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockSession
import nextflow.agent.rpc.AgentRpcRegistration

/**
 * White-box resume tests (design §7.3/§9.6, plan T4d/T7a/T7b/T7c/T7d): a cache hit
 * replays the stored generation through {@code collectOutputsV2} WITHOUT calling the
 * LLM. Built via the M1 task path ({@link AgentDef#buildAgentTask}) with a counting
 * mock runner; the replay is driven directly through {@code checkCachedOutput}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(30)
class AgentResumeIntegrationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        AgentCallInfo.clear()
    }

    /**
     * The minimal configuration that CONTAINERIZES a canonical agent task -- an enabled container
     * engine plus an `agent.container` image -- which every launch-spec runner now requires, since
     * its launch command is built from paths that exist only inside the runner image. An image
     * already present in the given config is preserved.
     */
    private static Map containerized(Map config = [:]) {
        final result = new LinkedHashMap(config)
        result.docker = [enabled: true]
        final agentScope = new LinkedHashMap((Map) (config.agent ?: [:]))
        agentScope.container = agentScope.container ?: 'agent-image:test'
        result.agent = agentScope
        return result
    }

    /** Spin a MockSession (MockExecutorFactory) and make it the global session. */
    private Session newSession(Map config = null) {
        def session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()
        Global.session = session
        return session
    }

    /** A probe observer counting cache-hit notifications. */
    private static TraceObserverV2 cachedProbe(List<TaskEvent> sink) {
        return new TraceObserverV2() {
            @Override void onTaskCached(TaskEvent event) { sink.add(event) }
        }
    }

    private static void injectObserver(Session session, TraceObserverV2 probe) {
        final f = Session.getDeclaredField('observersV2')
        f.setAccessible(true)
        final list = new ArrayList((List) f.get(session))
        list.add(probe)
        f.set(session, list)
    }

    private AgentDef newAgent(Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    /**
     * An agent owned by a MODULE script: {@code ScriptMeta.register} + {@code setScriptPath} is
     * what gives the owner a {@code moduleDir}, which is what makes {@code AgentDef.ownerBaseDir}
     * prefer it over {@code session.baseDir}. Without this the "module agent" and the plain local
     * agent would be the same object and these tests would prove nothing.
     */
    private AgentDef newModuleAgent(Path modDir, Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        final meta = ScriptMeta.register(owner)
        meta.setScriptPath(modDir.resolve('main.nf'))
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    /** Write `<modDir>/skills/<name>/SKILL.md` with the front matter SkillResolver expects. */
    private static Path moduleSkill(Path modDir, String name, String body) {
        final dir = modDir.resolve("skills/${name}")
        Files.createDirectories(dir)
        dir.resolve('SKILL.md').text = """\
            ---
            name: ${name}
            description: The ${name} skill
            ---

            ${body}
            """.stripIndent()
        return dir
    }

    private static AgentRunner namedRunner(String name = 'test') {
        return new AgentRunner() {
            @Override String getName() { name }
            @Override String run(AgentRunnerRequest req) { 'x' }
        }
    }

    /**
     * A legacy runner that records the requests it is handed. An anonymous class rather than a
     * coerced closure: {@code buildAgentTask} calls {@code getName()} on the selected runner, and a
     * closure proxy answers every interface method with the same closure body.
     */
    private static AgentRunner capturingRunner(List sink) {
        return new AgentRunner() {
            @Override String getName() { 'test' }
            @Override String run(AgentRunnerRequest req) { sink.add(req); 'ANSWER' }
        }
    }

    private static AgentRunner canonicalRunner(String name = 'canonical-test') {
        return new AgentRunner() {
            @Override String getName() { name }
            @Override AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(['/agent-rpc'], ['node'])
            }
            // the broker lives in the runner plugin, so a launch-spec runner must issue its own
            // registration; a canned one keeps this test off the network (see AgentRunner.register)
            @Override AgentRpcRegistration register(AgentRunnerRequest req, boolean remote) {
                new AgentRpcRegistration('inv-test', 'tok-test', '127.0.0.1:9999', 'fp-test')
            }
            @Override String run(AgentRunnerRequest req) { 'x' }
        }
    }

    def 'canonical agents default to local independently of process executor'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession(containerized([process: [executor: 'k8s']]))

        when:
        def processor = newAgent().buildAgentTask(['hello'])

        then:
        processor.config.executor == 'local'
    }

    def 'canonical agents resolve execution settings exclusively from agent scope'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        // `rpc.remoteHost` is required with a remote executor: the container is launched off the
        // driver host, so no container-engine host alias can stand in for the driver's address
        newSession(containerized([
            agent: [executor: 'k8s', container: 'agent-image:1', arch: 'arm64', cpus: 2, memory: '1 GB',
                    rpc: [remoteHost: 'driver.internal']],
            process: [executor: 'local', container: 'process-image:1', arch: 'amd64', cpus: 8, 'withName:qa': [executor: 'local']]
        ]))

        when:
        def processor = newAgent().buildAgentTask(['hello'])

        then:
        processor.config.executor == 'k8s'
        processor.config.container == 'agent-image:1'
        processor.config.arch == 'arm64'
        processor.config.cpus == 2
        processor.config.memory.toString() == '1 GB'
    }

    def 'canonical agent body resolves launch helpers statically under DELEGATE_ONLY'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession(containerized())
        def processor = newAgent().buildAgentTask(['hello'])
        processor.createStateObj()
        def task = processor.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        def body = ((Closure) processor.getTaskBody().closure.clone())
        body.setDelegate(new TaskContext(processor, [q: 'hello', task: task.config]))

        when:
        def command = body.call()

        then:
        command.startsWith("exec '/agent-rpc' '--endpoint'")
        command.contains("'--' 'node'")
    }

    def 'legacy runners reject a remote agent executor'() {
        given:
        AgentRunnerProvider.testRunner = namedRunner('legacy')
        newSession([agent: [executor: 'k8s']])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        def error = thrown(nextflow.exception.ScriptRuntimeException)
        error.message.contains('does not support executor `k8s`')
    }

    def 'legacy runners reject a SELECTOR-provided remote agent executor'() {
        given: 'the guard must read the RESOLVED executor, not just the plain scope'
        AgentRunnerProvider.testRunner = namedRunner('legacy')
        newSession([agent: ['withName:qa': [executor: 'k8s']]])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        def error = thrown(nextflow.exception.ScriptRuntimeException)
        error.message.contains('does not support executor `k8s`')
    }

    def 'an agent selector can opt a tool-free agent out of resume'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession(containerized([agent: ['withName:qa': [cache: false]]]))

        when:
        def processor = newAgent().buildAgentTask(['hello'])

        then:
        processor.config.isCacheable() == false
    }

    /**
     * Drive a cache-hit replay: build the (already-created) processor, make a TaskRun,
     * feed a TaskEntry whose context holds the output under `outName`, and assert the
     * replay binds the stored value with zero runner calls.
     */
    private Map replay(TaskProcessor processor, String outName, Object storedValue, AtomicInteger runnerCalls, List<TaskEvent> cachedEvents) {
        processor.createStateObj()
        final task = processor.createTaskRun(new TaskStartParams(TaskId.of(1), 1))

        final ctx = new TaskContext(processor, [(outName): storedValue])
        final entry = new TaskEntry(new TraceRecord(), ctx)
        final folder = Files.createTempDirectory('nxf-resume')
        final hash = CacheHelper.hasher('agent-resume').hash()

        final hit = processor.checkCachedOutput(task, folder, hash, entry)
        final bound = processor.getConfig().getOutputs().getParams()[0].getChannel().val
        return [hit: hit, bound: bound, runnerCalls: runnerCalls.get(), cachedEvents: cachedEvents.size(), task: task]
    }

    // -- T7a: a cache hit replays the stored generation through collectOutputsV2,
    //         binds the stored value, and makes ZERO runner (LLM) calls.
    def 'agent cache hit replays the stored output without calling the LLM'() {
        given:
        def runnerCalls = new AtomicInteger()
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> runnerCalls.incrementAndGet(); 'FRESH' } as AgentRunner
        def session = newSession()
        def cached = []
        injectObserver(session, cachedProbe(cached))

        and: 'the agent processor built via the M1 task path'
        def processor = newAgent().buildAgentTask(['hello'])

        when:
        def r = replay(processor, 'answer', 'STORED', runnerCalls, cached)

        then: 'checkCachedOutput reports a hit'
        r.hit == true
        and: 'the output channel replays the STORED value (collectOutputsV2 routed on the hit)'
        r.bound == 'STORED'
        and: 'the LLM runner was NEVER called (pure memoization, not a fresh completion)'
        r.runnerCalls == 0
        and: 'the cache-hit notification fired'
        r.cachedEvents == 1
        r.task.cached == true
    }

    // -- T7a (read guard, load-bearing): a cached entry with a NULL context must force a
    //         cache MISS. This makes the checkCachedOutput read guard `hasCacheableValues()
    //         && !entry.context` depend on the exec-type hasCacheableValues() branch
    //         (neutralizing that branch would let this falsely "hit").
    def 'a null cached context forces a cache miss (read guard depends on hasCacheableValues)'() {
        given:
        def runnerCalls = new AtomicInteger()
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> runnerCalls.incrementAndGet(); 'FRESH' } as AgentRunner
        newSession()

        and: 'the agent processor built via the M1 task path'
        def processor = newAgent().buildAgentTask(['hello'])
        processor.createStateObj()
        def task = processor.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        def folder = Files.createTempDirectory('nxf-resume-miss')
        def hash = CacheHelper.hasher('agent-resume-miss').hash()

        when: 'the cached entry has a NULL context (e.g. a pre-fix / cache-false entry)'
        def hit = processor.checkCachedOutput(task, folder, hash, new TaskEntry(new TraceRecord(), null))

        then: 'the exec body makes hasCacheableValues() true, so a missing context forces a miss'
        task.hasCacheableValues() == true
        hit == false
    }

    // -- T7a (write side): a real agent body run WRITES the declared output into the task
    //         context AND the exec-type hasCacheableValues() branch is true, so CacheDB.writeTaskEntry0's
    //         `proc.isCacheable() && task.hasCacheableValues()` gate would persist a non-null context.
    def 'agent body writes a non-null context and the task is cacheable (write-side gate)'() {
        given:
        def runnerCalls = new AtomicInteger()
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> runnerCalls.incrementAndGet(); 'ANSWER' } as AgentRunner
        newSession()

        and: 'the agent processor built via the M1 task path'
        def processor = newAgent().buildAgentTask(['hello'])
        processor.createStateObj()
        def task = processor.createTaskRun(new TaskStartParams(TaskId.of(1), 1))

        when: 'the synthetic GROOVY body runs against a task context seeded with the input'
        def body = ((Closure) processor.getTaskBody().closure.clone())
        def ctx = new TaskContext(processor, [q: 'hello'])
        body.setDelegate(ctx)
        body.call()

        then: 'the declared output landed in the context (CacheDB would persist it) and the runner ran once'
        ctx.get('answer') == 'ANSWER'
        runnerCalls.get() == 1
        and: 'the write-side gate operand hasCacheableValues() is true for the exec agent body'
        task.hasCacheableValues() == true
        processor.getConfig().isCacheable() == true
    }

    // -- T7c: `agent.cache = false` makes the agent non-cacheable, so no
    //         context is stored/looked up (the resume opt-out, design §7.5/D7).
    def 'cache false makes the agent non-cacheable (resume opt-out inherited)'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'x' } as AgentRunner
        newSession([agent: [cache: false]])

        when:
        def processor = newAgent().buildAgentTask(['hello'])

        then: 'the applied `cache false` directive disables caching -> no context storage/lookup'
        processor.getConfig().isCacheable() == false
    }

    // -- T7d: plain-`val` V2 process (NO agent) resumes correctly through the same shared
    //         cache path — bounds the blast radius of the hasCacheableValues() change.
    def 'plain val V2 process resumes (replays the stored value) - blast-radius regression'() {
        given:
        def session = newSession()
        def cached = []
        injectObserver(session, cachedProbe(cached))
        def owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }

        and: 'a plain V2 process with a single val output (no agent involved)'
        def config = new ProcessConfigV2(owner, 'plain')
        config.getOutputs().addParam('x', String, { getProperty('x') })
        config.getOutputs().getParams().each { it.setChannel(nextflow.extension.CH.value()) }
        def body = new BodyDef({ -> getDelegate().put('x', 'hello'); return null }, 'plain-src', 'exec')
        def processor = ProcessDef.createTaskProcessor(session, owner, 'plain', 'plain', 'plain', config, body)

        when:
        def r = replay(processor, 'x', 'hello', new AtomicInteger(), cached)

        then:
        r.hit == true
        r.bound == 'hello'
        r.cachedEvents == 1
    }

    // -- T4d: the synthetic BodyDef carries the canonical source + the prompt's valRefs.
    def 'synthetic BodyDef carries canonical source and folds prompt valRefs'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()
        def owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        def refs = [new nextflow.script.TokenValRef('params.threshold')]
        def agent = new AgentDef(owner, 'qa', [model: 'openai/gpt-4o'] as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q', refs))

        when:
        def processor = agent.buildAgentTask(['hello'])
        def body = processor.getTaskBody()

        then: 'source equals the canonical identity string built from effective values (single scalar output => null schema, default maxIter=20)'
        body.source == agent.canonicalAgentSource('openai/gpt-4o', 20, null, null, 'test')
        body.source.startsWith('agentRunner=test\nagentModel=openai/gpt-4o')
        and: 'the prompt free-var refs are folded into the BodyDef (so params.* enters the hash)'
        (body.getValNames() as Set) == ['params.threshold'] as Set
    }

    // -- T7b: two tasks differing ONLY in BodyDef.source produce different task hashes
    //         (proxy for "changing agent.model invalidates the cache").
    def 'changing the canonical source changes the task hash'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()
        def processor = newAgent().buildAgentTask(['hello'])
        processor.createStateObj()

        when:
        def t1 = processor.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        t1.source = 'agentModel=openai/gpt-4o\ntemperature=default'
        def t2 = processor.createTaskRun(new TaskStartParams(TaskId.of(2), 2))
        t2.source = 'agentModel=openai/gpt-4o-mini\ntemperature=default'
        def h1 = new nextflow.processor.TaskHasher(t1).compute()
        def h2 = new nextflow.processor.TaskHasher(t2).compute()

        then:
        h1 != h2
    }

    // -- ENDPOINT/CREDENTIAL IDENTITY (design D5). The RESOLVED endpoint is part of what an agent
    //    IS (a different endpoint serves a different model under the same id); the credential is
    //    not, so rotating a key must not invalidate a single stored entry -- nor be written into
    //    the key, which is a hash input persisted in the cache db.
    def 'the resolved endpoint enters the cache key; a rotated credential does not'() {
        given:
        AgentRunnerProvider.testRunner = namedRunner()
        and: 'an empty environment, so an exported OPENAI_* cannot supply a tier of its own'
        SysEnv.push([:])

        when: 'neither the endpoint nor the credential resolves'
        newSession()
        def plain = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'the key is the historical string -- no `baseUrl=` line, so existing entries stay valid'
        !plain.contains('baseUrl=')

        when: 'only the credential is configured, then rotated'
        newSession([agent: [apiKey: 'sk-first-2b8e']])
        def k1 = newAgent().buildAgentTask(['hello']).getTaskBody().source
        newSession([agent: [apiKey: 'sk-second-77c1']])
        def k2 = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'the key is unchanged by either value, and neither value appears in it'
        k1 == plain
        k2 == plain
        !k1.contains('sk-first-2b8e')

        when: 'an endpoint is configured'
        newSession([agent: [baseUrl: 'http://localhost:8000/v1']])
        def local = newAgent().buildAgentTask(['hello']).getTaskBody().source
        newSession([agent: [baseUrl: 'http://localhost:9000/v1']])
        def other = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'the key changes, so -resume re-executes rather than replaying another endpoint answer'
        local != plain
        local == plain + '\nbaseUrl=http://localhost:8000/v1'
        other != local

        when: 'the endpoint comes from the environment instead of the config'
        SysEnv.push([NXF_AGENT_BASE_URL: 'http://localhost:8000/v1'])
        newSession()
        def fromEnv = newAgent().buildAgentTask(['hello']).getTaskBody().source
        SysEnv.pop()

        then: 'the RESOLVED value is folded in, so the env tier keys identically to the config tier'
        fromEnv == local

        cleanup:
        SysEnv.pop()
    }

    // -- PROVIDER IDENTITY (design D6). An EXPLICIT `agent.apiProvider` selects which environment
    //    variables the endpoint and the credential come from, so it is part of how this agent was
    //    configured. An INFERRED one is a pure function of `baseUrl`, which is already in the key.
    def 'an explicit apiProvider enters the cache key; an inferred one does not'() {
        given:
        AgentRunnerProvider.testRunner = namedRunner()
        and: 'an empty environment, so an exported provider variable cannot supply a tier of its own'
        SysEnv.push([:])

        when: 'neither option is set'
        newSession()
        def plain = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'the historical string -- no `apiProvider=` line, so existing entries stay valid'
        !plain.contains('apiProvider=')

        when: 'the namespace is INFERRED from a well-known endpoint host'
        newSession([agent: [baseUrl: 'https://openrouter.ai/api/v1']])
        def inferred = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'only the endpoint enters: the inference adds no information the key does not have,'
        // and leaving it out means a later addition to the D3 host table cannot silently invalidate
        // a stored run of a pipeline nobody touched
        inferred == plain + '\nbaseUrl=https://openrouter.ai/api/v1'
        !inferred.contains('apiProvider=')

        when: 'the same endpoint is written together with an EXPLICIT namespace'
        newSession([agent: [baseUrl: 'https://openrouter.ai/api/v1', apiProvider: 'openrouter']])
        def explicit = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'the key changes: which variables were read IS part of how this agent was configured'
        explicit == inferred + '\napiProvider=openrouter'

        when: 'the namespace alone is set, redundantly to what the model prefix already implied'
        newSession([agent: [apiProvider: 'openai']])
        def redundant = newAgent().buildAgentTask(['hello']).getTaskBody().source

        then: 'still the documented one-time invalidation -- the key records the CONFIG, not the resolution'
        redundant == plain + '\napiProvider=openai'

        cleanup:
        SysEnv.pop()
    }

    def 'the widened ladder hands a provider-tier credential to the runner and persists none of it'() {
        given: 'an exported ANTHROPIC_API_KEY now resolves in core where the openai carve-out never'
        // reached it (design D2). It must reach the runner as a value and nothing else.
        def seen = []
        AgentRunnerProvider.testRunner = capturingRunner(seen)
        final secret = 'sk-ant-canary-6b12'
        SysEnv.push([ANTHROPIC_API_KEY: secret])
        newSession()

        when:
        def processor = newAgent([model: 'anthropic/claude-sonnet-4']).buildAgentTask(['hello'])
        def body = ((Closure) processor.getTaskBody().closure.clone())
        body.setDelegate(new TaskContext(processor, [q: 'hello']))
        body.call()

        then: 'the runner is handed the credential, already scoped to the model\'s own provider'
        seen.size() == 1
        seen[0].apiKey == secret

        and: 'while the cache key holds nothing of it -- rotating the key must replay, not re-run'
        !processor.getTaskBody().source.contains(secret)
        and: 'nor does the lineage record, nor the config the Platform observer serializes'
        !processor.config.get(AgentTaskInfo.CONFIG_KEY).toString().contains(secret)
        !processor.config.toString().contains(secret)

        cleanup:
        SysEnv.pop()
    }

    def 'the resolved credential reaches no persisted agent artifact'() {
        given:
        AgentRunnerProvider.testRunner = namedRunner()
        SysEnv.push([:])
        and:
        final secret = 'sk-leak-canary-3d91'
        newSession([agent: [apiKey: secret, baseUrl: 'http://localhost:8000/v1']])

        when:
        def processor = newAgent().buildAgentTask(['hello'])

        then: 'the resolved identity attached for the lineage observer holds no credential'
        def info = processor.config.get(AgentTaskInfo.CONFIG_KEY)
        info instanceof AgentTaskInfo
        !info.toString().contains(secret)

        and: 'nor does BodyDef.source, which TaskHasher folds into the hash stored in the cache db'
        !processor.getTaskBody().source.contains(secret)
        and: 'while the endpoint IS there -- it is not a secret and it IS part of the identity'
        processor.getTaskBody().source.contains('baseUrl=http://localhost:8000/v1')

        and: 'nor the resolved task config, which the Platform observer serializes as configText'
        !processor.config.toString().contains(secret)

        cleanup:
        SysEnv.pop()
    }

    // -- MODULE AGENT RESUME (agent-module design §4.5). canonicalAgentSource folds a name-sorted
    //    content fingerprint of every skill but NO agent name and NO file path, so a module agent's
    //    cache key tracks what its skills SAY, not where the module lives.

    def 'editing a module skill invalidates the agent cache key, moving the module does not'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()
        final root = Files.createTempDirectory('test')

        and: 'two module dirs holding a byte-identical skill'
        final modA = root.resolve('a'); Files.createDirectories(modA)
        moduleSkill(modA, 'greeting', 'Always greet the user by name.')
        final modB = root.resolve('b'); Files.createDirectories(modB)
        moduleSkill(modB, 'greeting', 'Always greet the user by name.')

        when:
        final srcA = newModuleAgent(modA, [model: 'openai/gpt-4o', skills: 'greeting'])
            .buildAgentTask(['hello']).getTaskBody().source
        final srcB = newModuleAgent(modB, [model: 'openai/gpt-4o', skills: 'greeting'])
            .buildAgentTask(['hello']).getTaskBody().source

        then: 'the skill fingerprint is in the key ...'
        srcA.contains('skills=')
        and: '... and a move/rename of the module dir does NOT invalidate it'
        srcA == srcB

        when: 'the module skill body is edited'
        moduleSkill(modB, 'greeting', 'Always greet the user by name, and sign off politely.')
        final srcB2 = newModuleAgent(modB, [model: 'openai/gpt-4o', skills: 'greeting'])
            .buildAgentTask(['hello']).getTaskBody().source

        then: 'the key changes, so -resume re-executes the agent'
        srcB2 != srcA
    }

    // -- MODULE-TOOL AGENT RESUME. A module/process-tool agent is NOT opted out of resume:
    //    canonicalAgentSource folds a fingerprint of every tool's
    //    descriptor AND backing script source, so a replay is only ever served for the exact tools
    //    that produced it. Without the fingerprint this would silently replay a stale generation
    //    after a tool edit -- which is why the blanket `cache false` existed.

    /** An in-scope typed process usable as an agent tool, owned by the agent's own script. */
    private AgentDef newToolAgent(String toolScript) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        final meta = ScriptMeta.register(owner)
        final config = new ProcessConfigV2(owner, 'uppercase')
        config.getInputs().addParam('text', String, false)
        config.getOutputs().addParam('result', String, { getProperty('result') })
        meta.addDefinition(new ProcessDef(owner, 'uppercase', config, new BodyDef({ -> toolScript }, toolScript, 'script')))
        return new AgentDef(owner, 'qa', [model: 'openai/gpt-4o', tools: 'nf:module_run:uppercase'] as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    def 'a module-tool agent is cacheable and its key tracks the tool script'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()

        when:
        def built = newToolAgent('tr a-z A-Z').buildAgentTask(['hello'])

        then: 'a module-tool agent is NOT opted out of resume ...'
        built.getConfig().isCacheable() == true
        and: '... because the tool identity is folded into the cache key'
        built.getTaskBody().source.contains('tools=')

        when: 'the tool process script is edited'
        def edited = newToolAgent('tr A-Z a-z').buildAgentTask(['hello']).getTaskBody().source

        then: 'the agent key changes, so -resume re-runs it instead of replaying a stale generation'
        edited != built.getTaskBody().source

        when: 'the tool is unchanged'
        def same = newToolAgent('tr a-z A-Z').buildAgentTask(['hello']).getTaskBody().source

        then: 'the key is stable, so -resume can replay'
        same == built.getTaskBody().source
    }

    def 'an fs:-tool agent participates in resume'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()

        when:
        def built = newAgent([model: 'openai/gpt-4o', tools: 'fs:*']).buildAgentTask(['hello'])

        then:
        built.getConfig().isCacheable() == true

        and: 'the capability is part of the key, so adding or dropping it re-runs the agent'
        def toolFree = newAgent([model: 'openai/gpt-4o']).buildAgentTask(['hello'])
        built.getTaskBody().source != toolFree.getTaskBody().source
    }

    def 'aliasing a module agent changes the task hash while leaving the canonical source alone'() {
        given:
        newSession()
        AgentRunnerProvider.testRunner = namedRunner()
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mod'); Files.createDirectories(mod)
        moduleSkill(mod, 'greeting', 'Always greet the user by name.')

        when:
        final declared = newModuleAgent(mod, [model: 'openai/gpt-4o', skills: 'greeting'])
        final aliased = (AgentDef) newModuleAgent(mod, [model: 'openai/gpt-4o', skills: 'greeting'])
            .cloneWithName('qc')
        final p1 = declared.buildAgentTask(['hello'])
        final p2 = aliased.buildAgentTask(['hello'])

        then: 'no agent name is folded into the canonical source'
        p1.getTaskBody().source == p2.getTaskBody().source
        and: 'but the alias renames the processor'
        p1.getName() == 'qa'
        p2.getName() == 'qc'

        when: 'the task hash folds the fully-qualified processor name'
        p1.createStateObj(); p2.createStateObj()
        final t1 = p1.createTaskRun(new TaskStartParams(TaskId.of(1), 1))
        t1.source = p1.getTaskBody().source
        final t2 = p2.createTaskRun(new TaskStartParams(TaskId.of(2), 2))
        t2.source = p2.getTaskBody().source

        then: 'aliasing therefore invalidates the resume cache'
        new nextflow.processor.TaskHasher(t1).compute() != new nextflow.processor.TaskHasher(t2).compute()
    }
}
