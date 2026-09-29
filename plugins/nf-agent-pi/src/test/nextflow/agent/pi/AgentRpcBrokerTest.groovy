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
package nextflow.agent.pi

import java.nio.charset.StandardCharsets
import java.security.MessageDigest
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ScheduledFuture
import java.util.concurrent.ScheduledThreadPoolExecutor
import java.util.concurrent.TimeUnit

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import io.grpc.Grpc
import io.grpc.ManagedChannel
import io.grpc.ManagedChannelBuilder
import io.grpc.MethodDescriptor
import io.grpc.ServerBuilder
import io.grpc.Status
import io.grpc.TlsChannelCredentials
import io.grpc.stub.ClientCalls
import io.grpc.stub.StreamObserver
import nextflow.Global
import nextflow.Session
import nextflow.agent.rpc.AgentRpcHost
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.AgentRunnerRequest
import nextflow.agent.ToolDescriptor
import nextflow.agent.ToolDispatcher
import org.slf4j.LoggerFactory
import spock.lang.Specification
import nextflow.agent.rpc.Probes

/**
 * Exercises the broker over a real gRPC connection. This is the guard for the
 * transport itself: the broker moved out of core into this plugin together with
 * its gRPC dependencies, so a resolution or packaging mistake shows up here as a
 * failure to start the server or to complete the handshake.
 */
class AgentRpcBrokerTest extends Specification {

    private static final MethodDescriptor.Marshaller<String> MARSHALLER = new MethodDescriptor.Marshaller<String>() {
        @Override InputStream stream(String value) {
            return new ByteArrayInputStream(value.getBytes(StandardCharsets.UTF_8))
        }
        @Override String parse(InputStream stream) {
            return new String(stream.readAllBytes(), StandardCharsets.UTF_8)
        }
    }

    // The full method name is spelled out as a LITERAL, not derived from
    // AgentRpcBroker.SERVICE_NAME/METHOD_NAME: agent-rpc/main.go hardcodes
    // `/nextflow.agent.AgentBroker/Connect`, so deriving it here would let a rename of
    // the constants keep this test green while breaking the Go proxy. The service name
    // deliberately keeps the `nextflow.agent` prefix even though the class now lives in
    // `nextflow.agent.pi` — it is a wire identifier, not a package reference.
    private static final MethodDescriptor<String,String> CONNECT = MethodDescriptor
        .<String,String>newBuilder()
        .setType(MethodDescriptor.MethodType.BIDI_STREAMING)
        .setFullMethodName('nextflow.agent.AgentBroker/Connect')
        .setRequestMarshaller(MARSHALLER)
        .setResponseMarshaller(MARSHALLER)
        .build()

    /**
     * The pre-connect budget the OLD rule produced: {@code max(requestTimeout,30) + 60} with the
     * default {@code agent.requestTimeout} of 120s. Any queueing executor -- k8s under pressure,
     * slurm, AWS Batch -- routinely waits longer than this before the job starts.
     */
    private static final long OLD_EXPIRY_SECONDS = 180

    AgentRpcBroker broker
    List<ManagedChannel> channels = []
    ListAppender<ILoggingEvent> logs
    Logger brokerLog

    def setup() {
        // attached before the first broker is built: the TLS opt-out warns from the constructor
        brokerLog = (Logger) LoggerFactory.getLogger(AgentRpcBroker)
        logs = new ListAppender<ILoggingEvent>()
        logs.start()
        brokerLog.addAppender(logs)
        driverSession([:])
        broker = AgentRpcBroker.get()
    }

    def cleanup() {
        channels.each { it.shutdownNow() }
        broker?.close()
        // a feature that lowered the level to observe the registration line must not leave it there
        brokerLog?.setLevel(null)
        brokerLog?.detachAppender(logs)
        Global.session = null
        AgentRpcHostResolver.reset()
    }

    /**
     * The session the broker resolves its advertised address against, with the DRIVER HOST placed
     * explicitly: an ordinary uncontainerized Linux host with one address. Without this the address
     * ladder would read the machine running the suite, and a suite running inside a container would
     * see {@code /.dockerenv}, take the containerized-driver row, and advertise that container's own
     * address where these specs expect the engine host alias.
     */
    private Session driverSession(Map config) {
        final session = new Session(config)
        AgentRpcHostResolver.install(session, new HostProbes())
        return Global.session = session
    }

    /** @see #driverSession */
    static class HostProbes implements Probes {
        @Override String outboundAddress() { '10.0.3.17' }
        @Override List<String> interfaceAddresses() { ['10.0.3.17'] }
    }

    /**
     * Replaces the singleton broker with one built from an explicit `agent.rpc` scope, so a test
     * can exercise a capability lifetime it can actually wait out.
     */
    private AgentRpcBroker rebuild(Map rpcOpts) {
        broker.close()
        driverSession([agent: [rpc: rpcOpts]])
        return broker = AgentRpcBroker.get()
    }

    def 'should complete the handshake and dispatch a brokered tool call'() {
        given:
        def dispatched = []
        def request = new AgentRunnerRequest(
            model: 'openai/test',
            prompt: 'say hello',
            maxIterations: 5,
            requestTimeoutSeconds: 30,
            agentName: 'demo',
            toolSpecs: [new ToolDescriptor('echo', 'echo it back', [type: 'object'], null)],
            dispatch: { String name, String args ->
                dispatched << [name, args]
                return '{"value":"hello"}'
            } as ToolDispatcher)
        def registration = broker.register(request, false)

        when:
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'the broker admits the capability and replies with the portable spec'
        start.type == 'start'
        start.protocolVersion == 2
        start.invocationId == registration.invocationId
        start.spec.model == 'openai/test'
        start.spec.prompt == 'say hello'

        when: 'an authorized tool is called'
        frames.send([
            type: 'tool_call',
            invocationId: registration.invocationId,
            callId: 'call-1',
            name: 'echo',
            arguments: [msg: 'hello'] ])
        def result = frames.next()

        then: 'the call reaches the dispatcher and the result comes back over the stream'
        result.type == 'tool_result'
        result.callId == 'call-1'
        result.isError == false
        result.result == '{"value":"hello"}'
        dispatched.size() == 1
        dispatched[0][0] == 'echo'
        dispatched[0][1].contains('hello')
    }

    def 'should refuse an unauthorized tool without reaching the dispatcher'() {
        given:
        def dispatched = []
        def request = new AgentRunnerRequest(
            model: 'openai/test',
            maxIterations: 5,
            requestTimeoutSeconds: 30,
            toolSpecs: [new ToolDescriptor('echo', 'echo it back', [type: 'object'], null)],
            dispatch: { String name, String args -> dispatched << name; return 'never' } as ToolDispatcher)
        def registration = broker.register(request, false)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        frames.next()

        when:
        frames.send([
            type: 'tool_call',
            invocationId: registration.invocationId,
            callId: 'call-1',
            name: 'not-declared',
            arguments: [:] ])
        def result = frames.next()

        then:
        result.type == 'tool_result'
        result.isError == true
        result.result.contains('not authorized')
        dispatched.isEmpty()
    }

    def 'should send the runner-native names in the start frame and authorize none of them'() {
        given: 'the two halves of the tool split: one brokered process, three runner-native leaves'
        def dispatched = []
        def request = new AgentRunnerRequest(
            model: 'openai/test',
            maxIterations: 5,
            requestTimeoutSeconds: 30,
            toolSpecs: [new ToolDescriptor('echo', 'echo it back', [type: 'object'], null)],
            nativeToolNames: ['read', 'write', 'bash'],
            dispatch: { String name, String args -> dispatched << name; return 'never' } as ToolDispatcher)
        def registration = broker.register(request, false)

        when:
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'the native names travel BESIDE the descriptors, so the harness can enable the'
        // matching pi builtins -- they are not descriptors and carry no schema of ours
        start.spec.nativeToolNames == ['read', 'write', 'bash']
        start.spec.toolSpecs*.name == ['echo']

        when: 'the model calls one of them over the RPC stream anyway'
        frames.send([
            type: 'tool_call',
            invocationId: registration.invocationId,
            callId: 'call-1',
            name: 'bash',
            arguments: [command: 'cat /etc/shadow'] ])
        def result = frames.next()

        then: 'the allowlist is built from the BROKERED half only, so a container-side tool cannot'
        // be relocated into the driver JVM by calling it back
        result.type == 'tool_result'
        result.isError == true
        result.result.contains('not authorized')
        dispatched.isEmpty()
    }

    def 'should refuse to register a request whose two tool halves overlap'() {
        given: 'a process named `read` alongside the `fs:read` the runner serves itself'
        def request = new AgentRunnerRequest(
            model: 'openai/test',
            requestTimeoutSeconds: 30,
            toolSpecs: [new ToolDescriptor('read', 'a process', [type: 'object'], null)],
            nativeToolNames: ['read'])

        when:
        broker.register(request, false)

        then: 'no capability is minted at all - the partition is checked before a job is submitted'
        def err = thrown(IllegalStateException)
        err.message.contains('partition violated')
    }

    def 'should reject a connection presenting a wrong token'() {
        given:
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)
        def frames = connect(registration.endpoint)

        when:
        frames.send([type: 'connect', invocationId: registration.invocationId, token: 'not-the-token'])

        then: 'the stream is failed instead of admitting the caller'
        frames.error() != null
    }

    def 'should keep an unconsumed capability valid far beyond the old requestTimeout-derived expiry'() {
        given: 'an invocation whose task sits queued -- register() runs when the SCRIPT is generated'
        def request = new AgentRunnerRequest(model: 'openai/test', prompt: 'queued for ages', requestTimeoutSeconds: 120)
        def registration = broker.register(request, false)

        expect: 'nothing can invalidate the capability for another hour, so a connect at t=180s+ finds it'
        // The delay is read off the armed deadline rather than waited out: the ONLY thing that removes
        // an unconsumed capability is that scheduled task, so "the job is still queued ten minutes
        // from now" is exactly "no expiry is due for an hour". Under max(requestTimeout,30)+60 this
        // same capability was dead at 180s and the task then failed with `Invalid agent RPC invocation
        // identity or token` -- a security-shaped message for what was only scheduler latency.
        pendingExpiry(registration.invocationId).getDelay(TimeUnit.SECONDS) > OLD_EXPIRY_SECONDS
        pendingExpiry(registration.invocationId).getDelay(TimeUnit.MINUTES) >= 59

        and: 'the per-request LLM timeout cannot shorten the capability either, at any value'
        pendingExpiry(broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 1), false)
            .invocationId).getDelay(TimeUnit.MINUTES) >= 59

        when: 'the executor finally releases the job and the proxy dials back'
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'the capability is still valid and the handshake completes'
        start.type == 'start'
        start.invocationId == registration.invocationId
        start.spec.prompt == 'queued for ages'
    }

    def 'should disarm the capability clock when the connect frame is accepted'() {
        given:
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 120), false)
        def expiry = pendingExpiry(registration.invocationId)

        expect: 'the deadline is armed at register(), because the PRE-connect wait is all it bounds'
        !expiry.isCancelled()
        !expiry.isDone()

        when:
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'consuming the capability cancels it, and nothing is armed in its place'
        // A live stream deliberately carries NO deadline. `requestTimeoutSeconds` bounds one LLM call,
        // so a legitimate agent -- maxIterations 20 x 120s, plus driver-side tool dispatch that itself
        // queues -- runs for tens of minutes; re-arming that constant here would kill nearly every
        // multi-turn agent mid-stream, which is the same spurious security-shaped failure D1 removes.
        // The next test proves a live stream outlives the capability budget.
        start.type == 'start'
        expiry.isCancelled()
    }

    def 'should dequeue a consumed capability deadline, not merely cancel it'() {
        given: 'two capabilities pending on the one-hour default budget'
        def first = broker.register(new AgentRunnerRequest(model: 'openai/test', prompt: 'one', requestTimeoutSeconds: 30), false)
        broker.register(new AgentRunnerRequest(model: 'openai/test', prompt: 'two', requestTimeoutSeconds: 30), false)

        expect:
        expiryQueueSize() == 2

        when: 'the first task starts and its proxy connects'
        def frames = connect(first.endpoint)
        frames.send([type: 'connect', invocationId: first.invocationId, token: first.token])
        frames.next()

        then: 'its deadline leaves the delay queue there and then'
        // ScheduledThreadPoolExecutor defaults removeOnCancelPolicy to false, so a cancelled task
        // sits in the queue until its ORIGINAL deadline -- and the queued Runnable captures the
        // Invocation, hence the prompt, the serialized inputs, the tool specs and the dispatcher
        // closure. Widening the pre-connect budget to an hour without this would mean an hour of
        // retention per *finished* agent, which is worst exactly where fan-out is widest.
        expiryQueueSize() == 1
    }

    def 'should honour a zero capability timeout rather than widening it to the default'() {
        given: 'the tightest window an operator can ask for'
        // nextflow.util.Duration is falsy at zero, so defaulting with Elvis here would silently
        // replace the most restrictive setting with the most permissive one -- and silently
        // relaxing a security knob is the one direction a default must never take.
        broker = rebuild([capabilityTimeout: '0s'])
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 3600), false)

        when:
        sleep(500)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])

        then: 'the capability was already gone, rather than valid for another hour'
        frames.error() != null
    }

    def 'should reject a connect frame once the capability timeout has elapsed'() {
        given: 'a capability that expires in a second, and a much longer per-request LLM timeout'
        broker = rebuild([capabilityTimeout: '1s'])
        // requestTimeoutSeconds is deliberately large: the pre-connect budget must come from
        // `agent.rpc.capabilityTimeout` ALONE. Under the previous max(requestTimeout,30)+60 rule
        // this capability would still be valid for another hour.
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 3600), false)

        when:
        sleep(1_500)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])

        then: 'the expired capability is gone from the broker and the stream is failed'
        frames.error() != null
    }

    def 'should not put a connected stream on the capability clock'() {
        given: 'a capability that expires in a second, consumed immediately'
        broker = rebuild([capabilityTimeout: '1s'])
        def request = new AgentRunnerRequest(
            model: 'openai/test',
            maxIterations: 5,
            requestTimeoutSeconds: 30,
            toolSpecs: [new ToolDescriptor('echo', 'echo it back', [type: 'object'], null)],
            dispatch: { String name, String args -> return '{"value":"hello"}' } as ToolDispatcher)
        def registration = broker.register(request, false)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        frames.next()

        when: 'the live stream outlives the capability timeout'
        sleep(1_500)
        frames.send([
            type: 'tool_call',
            invocationId: registration.invocationId,
            callId: 'call-1',
            name: 'echo',
            arguments: [:] ])
        def result = frames.next()

        then: 'the timeout bounds the PRE-connect wait only: an agent may legitimately run for hours'
        result.type == 'tool_result'
        result.isError == false
    }

    def 'should reject a replayed connect reusing a consumed capability'() {
        given:
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)
        def first = connect(registration.endpoint)
        first.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        first.next()

        when: 'a second stream presents the same identity and token'
        def replay = connect(registration.endpoint)
        replay.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])

        then: 'single-use consumption still rejects it, independently of any expiry'
        replay.error() != null
    }

    def 'should report the capabilities a resumed run leaves behind, on the default budget'() {
        given: 'fifty agent tasks whose scripts were generated and which then hit the resume cache'
        // TaskProcessor calls task.resolve(taskBody) -- hence register() -- BEFORE checkStoredOutput
        // and before checkCachedOrLaunchTask, so on `nextflow run -resume` a cache HIT still mints a
        // capability. Nothing connects with it and nothing cancels it.
        50.times {
            broker.register(new AgentRunnerRequest(model: 'openai/test', prompt: "cached ${it}", requestTimeoutSeconds: 30), false)
        }

        expect: 'no deadline has fired, and none is anywhere near due'
        // This is why counting EXPIRIES cannot surface the resume shape: on the default one-hour
        // budget such a run is over in seconds, so the fired-deadline count is zero and the operator
        // gets silence. The count that names it is what is still PENDING when the broker shuts down.
        pendingCount() == 50
        warnings().isEmpty()

        when: 'the run ends'
        broker.close()

        then: 'the retention is named and counted, rather than being released without a word'
        warnings().any {
            it.contains('50 agent capabilities registered but never consumed') && it.contains('resume cache')
        }

        and: 'reported once, though close() is reachable twice for the same broker'
        // cleanup() closes it again; session.onShutdown plus an explicit close is the same shape
        broker.close()
        warnings().count { it.contains('registered but never consumed') } == 1
    }

    def 'should tell an expired capability apart from an identity it never issued'() {
        given: 'a capability whose task takes longer to start than the budget allows'
        broker = rebuild([capabilityTimeout: '1s'])
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 3600), false)

        when: 'the proxy finally dials back'
        sleep(1_500)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def expired = frames.error()

        then: 'the task is told the budget lapsed, not that it presented a forged credential'
        // The reason has to travel as a real gRPC status: grpc-java closes a call aborted by a thrown
        // exception with `UNKNOWN: Application error processing RPC` and does not serialize the cause,
        // so every message this replaces was invisible at the far end.
        Status.fromThrowable(expired).code == Status.Code.UNAUTHENTICATED
        Status.fromThrowable(expired).description.contains('expired')
        Status.fromThrowable(expired).description.contains('agent.rpc.capabilityTimeout')

        and: 'the deadline that fired left evidence, naming the invocation and the budget'
        // it used to fire silently -- `{ invocations.remove(id, invocation) }` and nothing else
        warnings().any { it.contains(registration.invocationId) && it.contains('expired after 1s') }

        when: 'an identity the broker never issued is presented'
        def unknown = connect(registration.endpoint)
        unknown.send([type: 'connect', invocationId: UUID.randomUUID().toString(), token: registration.token])
        def rejected = unknown.error()

        then: 'there is nothing to disclose, so this one keeps the indiscriminate answer'
        Status.fromThrowable(rejected).code == Status.Code.UNAUTHENTICATED
        Status.fromThrowable(rejected).description == 'Invalid agent RPC invocation identity or token'
    }

    def 'should name the retry cause when an already-consumed capability comes back'() {
        given:
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)
        def first = connect(registration.endpoint)
        first.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        first.next()

        when: 'the same identity and token are presented a second time'
        def again = connect(registration.endpoint)
        again.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def error = again.error()

        then: 'the answer comes from the outcome record, because the entry itself is already gone'
        // The first connect REMOVED the entry, so this lookup finds null -- exactly what an id the
        // broker never issued finds. Without a record of how the capability left, the two are
        // indistinguishable and this path reports `Invalid ... identity or token`.
        Status.fromThrowable(error).code == Status.Code.UNAUTHENTICATED
        Status.fromThrowable(error).description.contains('already consumed')

        and: 'and it names the operational cause, not only the replay reading'
        // TaskProcessor re-submits a ProcessRetryableException / CloudSpotTerminationException copy
        // through task.makeCopy() WITHOUT re-resolving the task body, so a node termination or a spot
        // reclaim re-runs the SAME --invocation/--token. A first attempt that got as far as
        // connecting therefore lands here legitimately, and calling that a replay attack is the
        // security-shaped mislabel this branch exists to remove.
        Status.fromThrowable(error).description.contains('retried')
        warnings().any { it.contains(registration.invocationId) && it.contains('retried') }
    }

    def 'should refuse a connect frame that carries no invocation identity'() {
        given: 'a live capability, so the broker is serving and the pending map is not empty'
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)
        def frames = connect(registration.endpoint)

        when: 'a connect frame omits invocationId altogether'
        frames.send([type: 'connect', token: registration.token])
        def error = frames.error()

        then: 'the caller is refused with a status rather than with an opaque close'
        // `invocations` is a ConcurrentHashMap, so the null key threw NullPointerException out of
        // onNext and grpc-java closed the call with `UNKNOWN: Application error processing RPC`,
        // discarding the cause -- the exact failure this branch exists to stop emitting, and
        // reachable before authentication by any peer that can reach the port.
        Status.fromThrowable(error).code == Status.Code.UNAUTHENTICATED
        Status.fromThrowable(error).description == 'Invalid agent RPC invocation identity or token'

        and: 'the live capability is untouched, so the real proxy can still spend it'
        pendingCount() == 1
    }

    def 'should still name the consumed cause after a mass lapse of other capabilities'() {
        given: 'a budget long enough to connect within, and short enough that later capabilities lapse'
        broker = rebuild([capabilityTimeout: '2s'])
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)
        def first = connect(registration.endpoint)
        first.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        assert first.next().type == 'start'

        when: 'the run goes on to lapse more capabilities than one 1024-row budget could hold'
        1100.times { broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false) }
        awaitDrained()

        and: 'the connected task is retried, presenting the same capability a second time'
        def again = connect(registration.endpoint)
        again.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def error = again.error()

        then: 'it still learns an earlier attempt had connected, instead of being told it forged the token'
        // This CONSUMED row is the OLDEST, so a shared 1024-row budget let the 1100 EXPIRED rows evict
        // precisely the one a retry needs. The message then fell back to `Invalid agent RPC invocation
        // identity or token`: the security-shaped mislabel the record exists to remove, reappearing at
        // exactly the run sizes where retries are most likely.
        Status.fromThrowable(error).code == Status.Code.UNAUTHENTICATED
        Status.fromThrowable(error).description.contains('already consumed')
    }

    def 'should keep a wide fan-out diagnosable without burying the run in lapse warnings'() {
        given: 'a fan-out far wider than the old record held, every capability lapsing at once'
        broker = rebuild([capabilityTimeout: '0s'])

        when:
        1124.times { broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false) }
        awaitDrained()

        then: 'every outcome is still on the record, so any of these tasks can still be told why'
        // 1124 is deliberately just past the 1024 rows this record used to hold, because that bound
        // silently traded away the diagnosis it exists to provide: the evicted rows are the OLDEST,
        // i.e. the tasks whose retries arrive latest. The cap still exists at OUTCOME_HISTORY and is
        // two orders of magnitude up; a suite cannot reach it without a 100k-task run, so what is
        // asserted here is the property that was actually broken.
        outcomeCount() == 1124

        and: 'and a wide fan-out cannot bury the run in lapse warnings'
        warnings().count { it.contains('expired after') } == 10

        when: 'the broker shuts down'
        broker.close()

        then: 'the aggregate still carries the true total, whatever the per-lapse cap swallowed'
        warnings().any { it.contains('1124 agent capabilities that lapsed') }
    }

    def 'should bound a half-open agent stream with server keepalive'() {
        given:
        def builder = Mock(ServerBuilder)

        when:
        AgentRpcBroker.applyKeepAlive(builder)

        then: 'the server pings a silent connection, so a vanished node is noticed in ~80s not ~2h'
        // Without this there is NOTHING watching a connected stream: the post-connect deadline was
        // deliberately removed, an OOM-killed pod or a reclaimed spot instance never sends a FIN, and
        // both grpc-java's server default and Linux `tcp_keepalive_time` are 2h.
        1 * builder.keepAliveTime(60, TimeUnit.SECONDS)
        1 * builder.keepAliveTimeout(20, TimeUnit.SECONDS)

        and: 'the enforcement floor is grpc-go\'s own clamp, not a number agreed with the proxy'
        // grpc-go raises keepalive.ClientParameters.Time to KeepaliveMinPingTime = 10s whatever the
        // proxy configures, so a floor of 10s cannot be violated by any build of agent-rpc. A floor
        // stricter than the client's interval makes the server answer healthy streams with GOAWAY --
        // pinning it to the client library's minimum is what removes that coordination failure.
        1 * builder.permitKeepAliveTime(10, TimeUnit.SECONDS)

        and: 'permitKeepAliveWithoutCalls is left OFF, deliberately'
        // The server asks for no client certificate and binds every interface, so ping enforcement is
        // the only thing that evicts a peer which completes the TLS handshake and never opens a
        // stream. Turning it on to cover the pre-stream window would be a fail-open for a window that
        // is sub-millisecond: grpc-go dials lazily on the first RPC, and KeepAliveEnforcer resets its
        // strike count on stream creation.
        0 * builder.permitKeepAliveWithoutCalls(_)

        and: 'and nothing else is imposed on the connection'
        // maxConnectionAge / maxConnectionIdle would tear down a live agent mid-stream, which is the
        // spurious kill the pre/post-connect split removed in the first place.
        0 * builder._
    }

    def 'should serve TLS with a per-run certificate and pin it by its DER digest'() {
        given: 'the default broker, i.e. transport security on'
        def request = new AgentRunnerRequest(model: 'openai/test', prompt: 'confidential', requestTimeoutSeconds: 30)
        def registration = broker.register(request, false)

        expect: 'the pin is the SHA-256 of the served certificate DER, lowercase hex, no separators'
        // This is the whole cross-language contract: agent-rpc hashes rawCerts[0], which IS this DER.
        // Hashing the PEM text (or its base64 body) would produce an equally well-formed digest that
        // never matches, and the failure would be indistinguishable from a real pinning rejection.
        registration.fingerprint ==~ /[0-9a-f]{64}/
        registration.fingerprint == MessageDigest.getInstance('SHA-256')
            .digest(servedCertificateDer(broker.certificatePem)).encodeHex().toString()

        and: 'the proxy is told to pin it, rather than being left to infer trust from a missing flag'
        registration.transportArgs() == ['--fingerprint', registration.fingerprint]

        and: 'one identity per run, not per invocation'
        broker.register(request, false).fingerprint == registration.fingerprint

        when: 'a client that trusts only that certificate connects'
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'the handshake completes, so the server really serves the pinned certificate'
        start.type == 'start'
        start.spec.prompt == 'confidential'

        when: 'a cleartext client dials the same port'
        def plaintext = ManagedChannelBuilder.forTarget(registration.endpoint).usePlaintext().build()
        channels << plaintext
        // no frame is sent: the RPC alone drives the connection, and a cleartext h2 preface against
        // a TLS listener fails the handshake, so the failure arrives without touching the stream
        def rejected = framesOn(plaintext)

        then: 'there is no h2c fallback, so no frame -- and no prompt -- can cross in cleartext'
        rejected.error() != null
    }

    def 'should refuse a client that trusts a certificate the broker does not serve'() {
        given: 'the driver identity a task pins, and an impostor identity of exactly the same shape'
        def registration = broker.register(
            new AgentRunnerRequest(model: 'openai/test', prompt: 'confidential', requestTimeoutSeconds: 30), false)
        def impostor = AgentRpcTlsCredentials.create()

        expect: 'a different certificate is a different pin -- subject and issuer are identical'
        // the digest assertion is what stops this test from passing vacuously if the broker ever
        // stopped serving TLS: a cleartext listener also refuses a TLS client
        registration.fingerprint ==~ /[0-9a-f]{64}/
        impostor.fingerprint != registration.fingerprint

        when: 'the client trusts only the impostor, as a proxy handed the wrong digest would'
        final parts = registration.endpoint.split(':')
        final credentials = TlsChannelCredentials.newBuilder()
            .trustManager(new ByteArrayInputStream(impostor.certificatePem.getBytes(StandardCharsets.US_ASCII)))
            .build()
        final channel = Grpc.newChannelBuilderForAddress(parts[0], parts[1] as int, credentials).build()
        channels << channel
        def frames = framesOn(channel)

        then: 'the handshake fails, so no token and no start frame -- hence no prompt -- crosses'
        // This is the JVM half of pinning: a grpc-java client trusts "this exact certificate", so a
        // mismatch surfaces as a handshake failure. The requirement that the PROXY report a clear
        // pinning error rather than a generic TLS message is asserted where the digest comparison
        // actually lives -- agent-rpc's TestPinnedFingerprintRejectsAnotherCertificate.
        frames.error() != null
    }

    def 'should serve cleartext and say so explicitly when tls is disabled'() {
        given: 'the secure default broker built in setup() warned about nothing'
        assert warnings().isEmpty()

        when:
        broker = rebuild([tls: false])
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)

        then: 'the opt-out is loud: a cleartext run must not be an unnoticed one'
        warnings().any { it.contains('agent.rpc.tls=false') && it.contains('cleartext') }

        and: 'no certificate, no pin'
        broker.certificatePem == null
        registration.fingerprint == null

        and: 'the opt-out is an explicit flag: an absent --fingerprint must never mean "unpinned"'
        registration.transportArgs() == ['--insecure']

        when: 'a cleartext client connects'
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])

        then: 'the escape hatch still works'
        frames.next().type == 'start'
    }

    // -----------------------------------------------------------------------
    // In-band credential delivery (design D4). The credential core resolved travels BESIDE the
    // portable spec, as a top-level start-frame field, and only when the link is TLS-protected.
    // -----------------------------------------------------------------------

    def 'should carry the resolved credential beside the spec, not inside it'() {
        given: 'the default broker, i.e. transport security on'
        def registration = broker.register(new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            prompt: 'p',
            requestTimeoutSeconds: 30,
            apiKey: 'sk-in-band-7d41',
            baseUrl: 'https://gw.corp/v1'), false)

        when:
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'the credential is a TOP-LEVEL frame field, which runner.mjs installs with setRuntimeApiKey'
        start.apiKey == 'sk-in-band-7d41'

        and: 'and NOT inside the spec -- that half is what a transport may relay verbatim, log or persist'
        !((Map) start.spec).containsKey('apiKey')
        !start.spec.toString().contains('sk-in-band-7d41')

        and: 'the endpoint still travels on the spec, where it always did: it is not a secret'
        start.spec.baseUrl == 'https://gw.corp/v1'

        and: 'and nothing was warned -- TLS is on, so there is nothing to withhold'
        warnings().isEmpty()
    }

    def 'should send the resolved key and NEVER the no-credential placeholder'() {
        given: 'core resolved nothing while an endpoint IS declared -- exactly the shape that makes'
        // AgentRunnerRequest.credential() substitute `nxf-no-credential`. A runtime key OWNS its
        // provider in pi (setRuntimeApiKey is consulted before pi's own store and before the ambient
        // environment), so sending the placeholder here would SHADOW a credential an `env`/`secret`
        // channel had already delivered to the container and turn a working run into a 401.
        def request = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            prompt: 'p',
            requestTimeoutSeconds: 30,
            apiKey: null,
            baseUrl: 'https://gateway.corp/v1')

        expect: 'the placeholder is what credential() would have produced for this very request'
        request.credential() == AgentRunnerRequest.PLACEHOLDER_API_KEY

        when:
        def registration = broker.register(request, false)
        def frames = connect(registration.endpoint)
        frames.send([type: 'connect', invocationId: registration.invocationId, token: registration.token])
        def start = frames.next()

        then: 'and the frame carries no credential at all, under any key'
        start.type == 'start'
        !((Map) start).containsKey('apiKey')
        !start.toString().contains(AgentRunnerRequest.PLACEHOLDER_API_KEY)
    }

    def 'should withhold the credential when tls is disabled, and say so once'() {
        given: 'the cleartext escape hatch, where the frame is readable and the driver is not'
        // authenticated to the task -- shipping a credential there is strictly worse than today
        broker = rebuild([tls: false])
        final secret = 'sk-must-not-cross-2e88'
        def newRequest = { String prompt -> new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', prompt: prompt, requestTimeoutSeconds: 30, apiKey: secret) }

        when: 'two agent tasks connect over it'
        def first = broker.register(newRequest('one'), false)
        def f1 = connect(first.endpoint)
        f1.send([type: 'connect', invocationId: first.invocationId, token: first.token])
        def s1 = f1.next()
        def second = broker.register(newRequest('two'), false)
        def f2 = connect(second.endpoint)
        f2.send([type: 'connect', invocationId: second.invocationId, token: second.token])
        def s2 = f2.next()

        then: 'the escape hatch keeps working, it just carries no secret'
        s1.type == 'start'
        s1.spec.prompt == 'one'
        !((Map) s1).containsKey('apiKey')
        !((Map) s2).containsKey('apiKey')
        !s1.toString().contains(secret)

        and: 'the run is told why, and how to deliver the credential out of band instead'
        def withheld = warnings().findAll { it.contains('credential resolved by the driver is NOT sent') }
        withheld.size() == 1
        withheld[0].contains('agent.rpc.tls')
        withheld[0].contains('agent.containerOptions')
        withheld[0].contains('secret')

        and: 'ONCE per broker, not once per task: the cause is the session setting, not the invocation'
        // a wide fan-out would otherwise restate one fact hundreds of times over the run's real
        // errors -- the same reasoning as the per-lapse warning cap
        warnings().every { !it.contains(secret) }
    }

    def 'should warn -- never fail -- when the driver withheld a resolved credential'() {
        given: 'AgentConfig resolved a provider key and the endpoint gate refused to send it. The'
        // langchain4j runner treats this as fatal: core is its ONLY credential source. pi is not --
        // it reads its own auth store and the provider variables present in the CONTAINER, which
        // the driver cannot see and which is exactly how a Kubernetes Secret delivers a key. A
        // driver-side error here would break a deployment that is working, so it is a WARN.
        def newRequest = { String prompt -> new AgentRunnerRequest(
            model: 'openai/gpt-4o',
            prompt: prompt,
            requestTimeoutSeconds: 30,
            apiKey: null,
            apiProvider: 'openai',
            baseUrl: 'https://gw.corp/v1',
            credentialWithheld: true) }

        when: 'two tasks of the same fan-out connect'
        def first = broker.register(newRequest('one'), false)
        def f1 = connect(first.endpoint)
        f1.send([type: 'connect', invocationId: first.invocationId, token: first.token])
        def s1 = f1.next()
        def second = broker.register(newRequest('two'), false)
        def f2 = connect(second.endpoint)
        f2.send([type: 'connect', invocationId: second.invocationId, token: second.token])
        def s2 = f2.next()

        then: 'the run is NOT aborted -- both agents start, and pi resolves for itself'
        s1.type == 'start'
        s2.type == 'start'
        and: 'with no credential on the frame, and above all no placeholder'
        !((Map) s1).containsKey('apiKey')
        !((Map) s2).containsKey('apiKey')
        !s1.toString().contains(AgentRunnerRequest.PLACEHOLDER_API_KEY)

        and: 'one warning names the model, the endpoint, the namespace and the two remedies'
        def withheld = warnings().findAll { it.contains('was withheld from') }
        withheld.size() == 1
        withheld[0].contains('openai/gpt-4o')
        withheld[0].contains('https://gw.corp/v1')
        withheld[0].contains('agent.apiKey')
        withheld[0].contains('agent.apiProvider')

        and: 'said ONCE per model id, not once per task in the fan-out'
        warnings().size() == 1
    }

    /**
     * The address is resolved PER AGENT DEFINITION and rides on the request, so a run that mixes a
     * local-docker agent with, say, a Kubernetes one must advertise each task ITS OWN address. Before
     * this the broker held a single field and handed the first definition's address to every task,
     * which is the plausible-but-unroutable failure the whole ladder exists to prevent -- and the
     * cost of it is silent, because the mis-addressed task holds its capability for the full
     * `agent.rpc.capabilityTimeout`.
     */
    def 'each registration advertises the address resolved for ITS OWN agent'() {
        given:
        def alias = AgentRpcHost.of('host.docker.internal', 'docker host alias')
        def inCluster = AgentRpcHost.of('10.42.0.9', 'inferred from default route (in-cluster driver)')

        when: 'the local-docker agent registers first'
        def local = broker.register(remoteRequest(alias), true)
        def pod = broker.register(remoteRequest(inCluster), true)

        then: 'neither address leaks into the other registration'
        local.endpoint.startsWith('host.docker.internal:')
        pod.endpoint.startsWith('10.42.0.9:')
    }

    def 'the registration line carries the address and the ladder row that produced it'() {
        given: 'the line is DEBUG -- the source label is a diagnostic, not an operator warning'
        brokerLog.setLevel(Level.DEBUG)

        when:
        broker.register(remoteRequest(AgentRpcHost.of('10.0.3.17', 'docker host alias')), true)
        broker.register(remoteRequest(AgentRpcHost.of('driver.internal', 'agent.rpc.remoteHost')), true)

        then: 'an inferred address is indistinguishable from a configured one without this label'
        debugMessages().any { it ==~ /.*advertising 10\.0\.3\.17:\d+ \(docker host alias\).*/ }
        debugMessages().any { it ==~ /.*advertising driver\.internal:\d+ \(agent\.rpc\.remoteHost\).*/ }
    }

    def 'an address the ladder is not certain about is warned about once per ADDRESS'() {
        given:
        def multiHomed = AgentRpcHost.of('10.1.2.3', 'inferred from default route',
                ['the driver host is multi-homed: the default route selected `10.1.2.3`, but `192.168.40.7` also exist - set `agent.rpc.remoteHost`'])

        when: 'two tasks of the SAME agent register'
        broker.register(remoteRequest(multiHomed), true)
        broker.register(remoteRequest(multiHomed), true)

        then: 'the warning is a property of the address, so a fan-out must not repeat it per task'
        def warned = warnings().findAll { it.contains('multi-homed') }
        warned.size() == 1
        warned[0] ==~ /.*advertising 10\.1\.2\.3:\d+ - .*/
        warned[0].contains('192.168.40.7')
    }

    def 'a registration that carries no address of its own falls back to the session-level one'() {
        when: 'a runner that registers without going through the pre-ignition guard'
        def registration = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), true)

        then: 'shipped behaviour: the local rows, keyed off the session\'s enabled engine'
        registration.endpoint.startsWith('host.docker.internal:')
    }

    private static AgentRunnerRequest remoteRequest(AgentRpcHost host) {
        return new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30, brokerHost: host)
    }

    /** The broker warnings emitted so far in this feature method. */
    private List<String> warnings() {
        return new ArrayList<ILoggingEvent>(logs.list).findAll { it.level == Level.WARN }*.formattedMessage
    }

    /** @see #warnings */
    private List<String> debugMessages() {
        return new ArrayList<ILoggingEvent>(logs.list).findAll { it.level == Level.DEBUG }*.formattedMessage
    }

    /**
     * The pending pre-connect expiry of a capability that has not been consumed yet.
     *
     * Reaching into the broker is deliberate. The D1 regression is "a capability survives a queueing
     * delay of many minutes", and the honest way to observe that without a suite that takes minutes is
     * to read the deadline it is holding: nothing else removes an unconsumed capability, so an expiry
     * that is not due for an hour IS the guarantee.
     */
    private ScheduledFuture<?> pendingExpiry(String invocationId) {
        final field = AgentRpcBroker.getDeclaredField('invocations')
        field.setAccessible(true)
        final invocation = ((Map)field.get(broker)).get(invocationId)
        assert invocation != null : "No pending capability for invocation ${invocationId}"
        return (ScheduledFuture<?>) invocation.expiry
    }

    /** How many capability deadlines the broker is still holding, cancelled ones included. */
    private int expiryQueueSize() {
        final field = AgentRpcBroker.getDeclaredField('expiryPool')
        field.setAccessible(true)
        return ((ScheduledThreadPoolExecutor) field.get(broker)).getQueue().size()
    }

    /** How many capabilities are registered and still waiting for a connection. */
    private int pendingCount() {
        final field = AgentRpcBroker.getDeclaredField('invocations')
        field.setAccessible(true)
        return ((Map) field.get(broker)).size()
    }

    /** How many terminal outcomes the broker remembers, i.e. the size of the bounded record. */
    private int outcomeCount() {
        final field = AgentRpcBroker.getDeclaredField('outcomes')
        field.setAccessible(true)
        return ((Map) field.get(broker)).size()
    }

    /**
     * Waits for every pending capability to leave the map, by polling rather than by sleeping for a
     * guessed interval: the expiries run on the broker's own single-threaded pool, so how long a
     * thousand of them take is a property of the machine, not of the test.
     */
    private void awaitDrained() {
        final deadline = System.currentTimeMillis() + 20_000
        while( pendingCount() > 0 && System.currentTimeMillis() < deadline )
            sleep(10)
        assert pendingCount() == 0 : "Capabilities were still pending after 20s: ${pendingCount()}"
    }

    def 'should advertise the driver host a containerized task can reach the broker on'() {
        given: 'a session whose only container engine is the one under test'
        def registration = brokerWith(config).register(
            new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), true)

        expect: 'the ADVERTISED host, with the bound port - the server itself binds every interface'
        registration.endpoint ==~ /\Q${host}\E:\d+/

        where: 'the two engines that name their container host, and an explicit override'
        config                                                                  || host
        [docker: [enabled: true]]                                               || 'host.docker.internal'
        [podman: [enabled: true]]                                               || 'host.containers.internal'
        [docker: [enabled: true], agent: [rpc: [remoteHost: 'driver.svc']]]     || 'driver.svc'
        [singularity: [enabled: true], agent: [rpc: [remoteHost: '127.0.0.1']]] || '127.0.0.1'
    }

    def 'should refuse to advertise an unresolvable driver host'() {
        given: 'a microVM created with no network at all, which no address can reach'
        // singularity used to stand here, and no longer does: it creates no network namespace, so
        // the ladder now answers it with the driver's own loopback (error row E2 is what is left)
        brokerWith([smolvm: [enabled: true, network: false]])

        // AgentDef rejects this configuration before the run starts, but the runner SPI is public,
        // so the broker must not hand out an endpoint that reads literally `null:<port>`
        when: 'a remote registration is asked for anyway'
        broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), true)

        then: 'the ladder row itself is raised -- it names the fact that decided and the remedy, both of which a generic "no address is configured" line would discard'
        def e = thrown(IllegalStateException)
        e.message.contains('smolvm.network')

        when: 'the same broker registers a driver-local task'
        def local = broker.register(new AgentRunnerRequest(model: 'openai/test', requestTimeoutSeconds: 30), false)

        then: 'loopback needs no configuration'
        local.endpoint ==~ /127\.0\.0\.1:\d+/
    }

    /**
     * Replace the broker created by {@code setup()} with one built from the given session config, so
     * a spec can pin how the advertised endpoint is resolved. The broker reads its config once, in
     * its constructor, hence the rebuild.
     */
    private AgentRpcBroker brokerWith(Map config) {
        broker.close()
        driverSession(config)
        broker = AgentRpcBroker.get()
        return broker
    }

    /** Minimal blocking client over the broker's bidi JSON stream. */
    private Frames connect(String endpoint) {
        final parts = endpoint.split(':')
        final channel = channelFor(parts[0], parts[1] as int)
        channels << channel
        return framesOn(channel)
    }

    private Frames framesOn(ManagedChannel channel) {
        final inbound = new LinkedBlockingQueue<Object>()
        final observer = ClientCalls.asyncBidiStreamingCall(
            channel.newCall(CONNECT, io.grpc.CallOptions.DEFAULT),
            new StreamObserver<String>() {
                @Override void onNext(String value) { inbound.put(new JsonSlurper().parseText(value)) }
                @Override void onError(Throwable error) { inbound.put(error) }
                @Override void onCompleted() { inbound.put('completed') }
            })
        return new Frames(observer, inbound)
    }

    /**
     * Dials the broker the way the Go proxy does: no CA vouches for the per-run certificate, so the
     * only trust root is that certificate itself. Falls back to cleartext when the broker was built
     * with `agent.rpc.tls = false`, which is what the proxy's `--insecure` branch does.
     */
    private ManagedChannel channelFor(String host, int port) {
        final pem = broker.certificatePem
        if( !pem )
            return ManagedChannelBuilder.forAddress(host, port).usePlaintext().build()
        final credentials = TlsChannelCredentials.newBuilder()
            .trustManager(new ByteArrayInputStream(pem.getBytes(StandardCharsets.US_ASCII)))
            .build()
        return Grpc.newChannelBuilderForAddress(host, port, credentials).build()
    }

    /** The leaf certificate DER, i.e. the bytes both ends must hash to agree on a fingerprint. */
    private static byte[] servedCertificateDer(String pem) {
        final body = pem
            .replace('-----BEGIN CERTIFICATE-----', '')
            .replace('-----END CERTIFICATE-----', '')
            .replaceAll('\\s', '')
        return Base64.decoder.decode(body)
    }

    private static class Frames {
        private final StreamObserver<String> outbound
        private final LinkedBlockingQueue<Object> inbound

        Frames(StreamObserver<String> outbound, LinkedBlockingQueue<Object> inbound) {
            this.outbound = outbound
            this.inbound = inbound
        }

        void send(Map message) { outbound.onNext(JsonOutput.toJson(message)) }

        Map next() {
            final frame = take()
            if( !(frame instanceof Map) )
                throw new AssertionError("Expected a broker frame but received: ${frame}" as Object)
            return (Map)frame
        }

        Throwable error() {
            final frame = take()
            return frame instanceof Throwable ? (Throwable)frame : null
        }

        private Object take() {
            final frame = inbound.poll(20, TimeUnit.SECONDS)
            if( frame == null )
                throw new AssertionError('Timed out waiting for a broker frame' as Object)
            return frame
        }
    }
}
