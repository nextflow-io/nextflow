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
import java.security.SecureRandom
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledFuture
import java.util.concurrent.ScheduledThreadPoolExecutor
import java.util.concurrent.ThreadFactory
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicBoolean
import java.util.concurrent.atomic.AtomicInteger

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.grpc.BindableService
import io.grpc.MethodDescriptor
import io.grpc.Server
import io.grpc.ServerBuilder
import io.grpc.ServerServiceDefinition
import io.grpc.Status
import io.grpc.stub.ServerCalls
import io.grpc.stub.StreamObserver
import nextflow.Global
import nextflow.Session
import nextflow.agent.AgentProtocolSpec
import nextflow.agent.rpc.AgentRpcConfig
import nextflow.agent.rpc.AgentRpcHost
import nextflow.agent.rpc.AgentRpcRegistration
import nextflow.agent.AgentRunnerRequest

/**
 * Embedded, driver-side broker for canonical agent tasks.
 *
 * The POC deliberately uses a JSON gRPC marshaller so the semantic JSONL
 * protocol and the Go proxy can evolve before committing to generated message
 * classes. The service is still a real bidirectional gRPC stream.
 *
 * The broker and its gRPC dependencies live in this plugin rather than in core so
 * that a distribution with no agent runner installed carries no RPC stack. Nothing
 * here is Pi-specific: it only reads runner-neutral {@link AgentRunnerRequest}
 * fields, and it is hosted here because {@code nf-agent-pi} is currently the only
 * runner with a launch spec, and it already builds the matching Go proxy. Extract
 * it to a shared {@code nf-agent-rpc} plugin when a second canonical runner needs
 * it, rather than introducing plugin-to-plugin coupling now.
 */
@Slf4j
@CompileStatic
class AgentRpcBroker {

    static final String SERVICE_NAME = 'nextflow.agent.AgentBroker'
    static final String METHOD_NAME = 'Connect'

    private static final MethodDescriptor.Marshaller<String> STRING_MARSHALLER = new MethodDescriptor.Marshaller<String>() {
        @Override InputStream stream(String value) {
            return new ByteArrayInputStream(value.getBytes(StandardCharsets.UTF_8))
        }

        @Override String parse(InputStream stream) {
            return new String(stream.readAllBytes(), StandardCharsets.UTF_8)
        }
    }

    private static final MethodDescriptor<String,String> CONNECT_METHOD = MethodDescriptor
        .<String,String>newBuilder()
        .setType(MethodDescriptor.MethodType.BIDI_STREAMING)
        .setFullMethodName(MethodDescriptor.generateFullMethodName(SERVICE_NAME, METHOD_NAME))
        .setRequestMarshaller(STRING_MARSHALLER)
        .setResponseMarshaller(STRING_MARSHALLER)
        .build()

    /**
     * The endpoint advertised to a non-remote agent task, i.e. a plain child process in the driver's
     * own network namespace, for which loopback is always correct. Core no longer registers such a
     * task -- a canonical agent task is always containerized -- but the runner SPI still allows it.
     * This is only the advertised host, never the bind address: the server binds every interface.
     */
    private static final String LOCAL_HOST = '127.0.0.1'

    /**
     * How often the server pings an otherwise silent connection, and how long the peer then has to
     * answer. Removing the post-connect deadline is right -- an agent may legitimately run for hours
     * -- but on its own it leaves NOTHING watching a stream whose peer is simply gone: an OOM-killed
     * pod or a reclaimed spot instance disappears without ever sending a FIN, so the {@link
     * Invocation} and its {@link ResponseSink} would be held until the operating system's TCP
     * keepalive noticed, which on Linux is {@code tcp_keepalive_time} = 2h. grpc-java's own server
     * default is the same 2h, so bounding this needs an explicit setting.
     *
     * <p>60s + 20s detects a vanished node inside ~80s at the cost of a few bytes a minute. The ping
     * is a transport frame answered by the peer's networking loop, not by the application, so it does
     * NOT mistake "the model has been thinking for ten minutes" for a dead peer -- and the server
     * pings even while the stream is idle, because {@code NettyServerHandler} builds its
     * {@code KeepAliveManager} with {@code keepAliveDuringTransportIdle = true}. A silent-but-open
     * agent stream is exactly the case this has to cover.
     */
    private static final long KEEPALIVE_TIME_SECONDS = 60
    private static final long KEEPALIVE_TIMEOUT_SECONDS = 20

    /**
     * grpc-go clamps {@code keepalive.ClientParameters.Time} up to {@code KeepaliveMinPingTime} =
     * 10s, so no build of the proxy can ping the driver more often than this however its interval is
     * configured. The server's enforcement floor is derived from that CLIENT-LIBRARY floor rather
     * than from whatever number {@code agent-rpc} currently passes: a server that permits pings more
     * rarely than the client sends them answers with GOAWAY and kills a perfectly healthy stream, and
     * pinning the floor to the library's own minimum makes the two sides impossible to get out of
     * step when either is edited later.
     */
    private static final long PERMIT_KEEPALIVE_TIME_SECONDS = 10

    /**
     * How many terminal outcomes are remembered so an arriving {@code connect} can be told what
     * happened to a capability that is no longer pending, and how many expiries are reported at WARN
     * before the rest fall to DEBUG. Both exist to stop a diagnosis gap without opening a second
     * retention hole, so both are capped; the aggregate logged by {@link #reportRetention} carries
     * the true totals regardless of either cap.
     *
     * <p>The record has to be sized against the RUN, not against a guess at how many rejections are
     * interesting, because the row a retry needs is the one written when its FIRST attempt connected
     * -- arbitrarily far back. A 1024-row budget got this exactly wrong: it is shared between
     * {@code CONSUMED} and {@code EXPIRED}, so a wide fan-out or a cache-hit resume that lapses
     * thousands of capabilities evicts the oldest rows first, which are precisely the connected
     * tasks a node termination is most likely to retry. The retry then fell back to
     * {@code Invalid agent RPC invocation identity or token} -- the security-shaped mislabel the
     * record was added to remove, reappearing at the run sizes where retries are most likely.
     *
     * <p>One row is an id string and an enum reference, ~150 bytes, so this ceiling is ~15 MB and is
     * reached only by a run with more than 100k agent tasks. That is a fraction of what the same run
     * spends on the pending capabilities themselves, each of which pins its whole request.
     */
    private static final int OUTCOME_HISTORY = 100_000
    private static final int LAPSE_WARN_LIMIT = 10

    private static AgentRpcBroker singleton

    private final Map<String,Invocation> invocations = new ConcurrentHashMap<>()
    private final Map<String,Outcome> outcomes = Collections.synchronizedMap(new OutcomeHistory())
    private final AtomicInteger lapsedCount = new AtomicInteger()
    private final AtomicBoolean closed = new AtomicBoolean()
    /** Whether the "credential withheld, TLS is off" warning has already been emitted for this run. */
    private final AtomicBoolean insecureCredentialWarned = new AtomicBoolean()
    /**
     * The advertised addresses already warned about. The warnings are a property of the ADDRESS, so
     * emitting them per registration would repeat the same line once per agent task -- but a run
     * with several agent definitions can legitimately advertise several addresses (see
     * {@link AgentRunnerRequest#brokerHost}), so the gate is per address rather than per run.
     */
    private final Set<String> hostWarned = ConcurrentHashMap.newKeySet()
    /** Model ids whose "credential withheld by the endpoint gate" warning has already been emitted. */
    private final Set<String> withheldCredentialWarned = ConcurrentHashMap.newKeySet()
    private final ExecutorService dispatchPool = Executors.newCachedThreadPool(daemonFactory('nf-agent-rpc-dispatch'))
    private final ScheduledThreadPoolExecutor expiryPool = expiryPool()
    private final SecureRandom random = new SecureRandom()
    private final Server server
    /**
     * The address advertised to a remote agent task whose request carries NONE of its own -- a
     * runner that registers without the pre-ignition guard, or a spec. It can be an error row, whose
     * message {@link #register} raises verbatim. @see AgentRunnerRequest#brokerHost
     */
    private final AgentRpcHost fallbackHost
    private final long capabilitySeconds
    /** SHA-256 of the served certificate, lowercase hex; {@code null} when TLS is disabled. */
    private final String fingerprint
    /** Whether the broker deliberately serves cleartext, i.e. {@code agent.rpc.tls = false}. */
    private final boolean insecureTransport
    /** The served certificate, for a test client that has to trust it; public material. */
    private final String certificatePem

    private static ThreadFactory daemonFactory(String name) {
        return { Runnable task ->
            final thread = new Thread(task, name)
            thread.daemon = true
            return thread
        } as ThreadFactory
    }

    /**
     * Built directly rather than via {@code Executors.newSingleThreadScheduledExecutor}, whose
     * wrapper hides {@code setRemoveOnCancelPolicy}. The default policy is {@code false}, so a
     * cancelled expiry stays in the delay queue until its original deadline -- and the queued
     * {@link Runnable} pins the whole {@link Invocation}: the prompt, the serialized inputs, the
     * tool specs and the dispatcher closure. With the pre-connect budget now an hour, every
     * successfully consumed invocation would hold that for an hour after it finished. Removing on
     * cancel is what makes consumption actually release the capability.
     */
    private static ScheduledThreadPoolExecutor expiryPool() {
        final pool = new ScheduledThreadPoolExecutor(1, daemonFactory('nf-agent-rpc-expiry'))
        pool.setRemoveOnCancelPolicy(true)
        return pool
    }

    /**
     * How a capability left the pending map. A {@code connect} that arrives once the entry is gone
     * finds only {@code null} and cannot otherwise tell an hour-old queueing delay from a forged
     * identity, so both answers used to be the same security-shaped
     * {@code Invalid agent RPC invocation identity or token} -- the exact message this PR set out to
     * stop emitting for a non-security cause.
     */
    private static enum Outcome { CONSUMED, EXPIRED }

    /**
     * Terminal outcomes, newest last, capped at {@link #OUTCOME_HISTORY}. Evicting the oldest row is
     * safe here in the way that evicting the oldest PENDING capability would not be: losing a row
     * only downgrades a diagnostic message back to the generic rejection, where evicting a pending
     * capability would invalidate a task that is about to start.
     *
     * <p>Keying it by invocation id makes it a status oracle for anyone who knows an id -- and the id
     * is NOT secret: it is on the task's argv as {@code --invocation}, hence in {@code .command.run},
     * in the trace, and in {@code ps} on the node. Knowing "expired" versus "consumed" without the
     * token buys such a reader nothing, while telling the operator which one it was is the whole
     * point of keeping the record.
     */
    private static class OutcomeHistory extends LinkedHashMap<String,Outcome> {
        @Override protected boolean removeEldestEntry(Map.Entry<String,Outcome> eldest) {
            return size() > OUTCOME_HISTORY
        }
    }

    private static class Invocation {
        String id
        String token
        AgentRunnerRequest request
        Set<String> callIds = ConcurrentHashMap.<String>newKeySet()
        AtomicInteger callCount = new AtomicInteger()
        Set<String> allowedTools = Collections.emptySet()
        /**
         * The pending pre-connect expiry task, cancelled once the capability is consumed so a
         * connected stream leaves nothing scheduled. Written by {@code register} on the task-body
         * thread and read by the gRPC thread that accepts the {@code connect} frame.
         */
        volatile ScheduledFuture<?> expiry
    }

    /**
     * Serializes writes to a single response stream and guards against writing
     * after the stream has terminated. gRPC {@link StreamObserver}s are not
     * thread-safe and reject {@code onNext} after a terminal event, so pool
     * threads must funnel every send through here.
     */
    private static class ResponseSink {
        private final StreamObserver<String> responses
        private boolean closed

        ResponseSink(StreamObserver<String> responses) {
            this.responses = responses
        }

        void send(Map message) {
            synchronized(responses) {
                if( closed )
                    return
                responses.onNext(JsonOutput.toJson(message))
            }
        }

        void complete() {
            synchronized(responses) {
                if( closed )
                    return
                closed = true
                responses.onCompleted()
            }
        }

        void markClosed() {
            synchronized(responses) {
                closed = true
            }
        }

        /**
         * Fails the call with an explicit gRPC status, which is the ONLY way a reason reaches the
         * task. Throwing out of {@code onNext} cannot carry one: grpc-java's
         * {@code ServerImpl$JumpToApplicationThreadServerStreamListener.internalClose} closes with
         * {@code Status.UNKNOWN.withDescription("Application error processing RPC").withCause(t)},
         * and the cause is never serialized -- so every message this replaces reached the proxy as a
         * bare {@code UNKNOWN} and nothing else.
         */
        void fail(Status status) {
            synchronized(responses) {
                if( closed )
                    return
                closed = true
                responses.onError(status.asRuntimeException())
            }
        }
    }

    private AgentRpcBroker(AgentRpcConfig config, Session session) {
        final requestedPort = config.port
        // The address a registration actually advertises rides ON THE REQUEST: AgentDef's
        // pre-ignition guard resolves it per agent definition with the full context -- the executor
        // instance, the engine config, the task's container options -- none of which is recoverable
        // here, and a run may hold several (a local-docker agent and a k8s agent resolve different
        // addresses, and each task must be told its own). What is resolved here is only the FALLBACK
        // for a registration that carries none, and it is deliberately weaker: session-level facts
        // answer `docker` even for a run whose agent task is a Kubernetes pod.
        this.fallbackHost = config.resolveBrokerHost(session)
        // Explicitly null-checked, not Elvis: nextflow.util.Duration is falsy at zero, so `?:` would
        // silently widen `agent.rpc.capabilityTimeout = '0s'` -- the tightest window an operator can
        // ask for -- into the one-hour default, the wrong direction for a security-relevant knob.
        this.capabilitySeconds = (config.capabilityTimeout != null ? config.capabilityTimeout : AgentRpcConfig.DEFAULT_CAPABILITY_TIMEOUT).seconds
        this.insecureTransport = !config.tlsEnabled()
        final ServerBuilder<?> builder = ServerBuilder.forPort(requestedPort)
        builder.addService(service())
        applyKeepAlive(builder)
        if( config.tlsEnabled() ) {
            // A per-run, self-signed identity the task pins by fingerprint: it closes the payload
            // exposure (the start frame carries the prompt and inputJson, then every tool argument
            // and result crosses the same link) and it authenticates the driver, so a process that
            // occupies the advertised endpoint cannot serve a forged start frame. Both PEMs are fed
            // in as in-memory streams -- no key material reaches disk.
            final credentials = AgentRpcTlsCredentials.create()
            builder.useTransportSecurity(credentials.certificateStream(), credentials.privateKeyStream())
            this.certificatePem = credentials.certificatePem
            this.fingerprint = credentials.fingerprint
        }
        else {
            this.certificatePem = null
            this.fingerprint = null
            log.warn "Agent RPC transport security is disabled (agent.rpc.tls=false) -- the agent prompt, inputs, tool arguments and results cross the network in cleartext, and the driver is not authenticated to the task"
        }
        this.server = builder.build().start()
        session?.onShutdown { close() }
        // The fingerprint is a public commitment, not a secret, and logging it is what lets an
        // operator compare it against the two digests the proxy prints on a pinning failure.
        log.debug "Agent RPC broker listening on port ${server.port}${fingerprint ? " tls-fingerprint=${fingerprint}" : ' (cleartext)'}"
    }

    /**
     * The only liveness bound on a CONNECTED stream. See {@link #KEEPALIVE_TIME_SECONDS} for why a
     * connected stream cannot be put back on a deadline and what the keepalive replaces.
     *
     * <p>{@code permitKeepAliveWithoutCalls} is deliberately LEFT at its default of {@code false}.
     * It would only matter in the window before the proxy's first stream opens, and that window is
     * sub-millisecond -- grpc-go dials lazily on the first RPC, and {@code KeepAliveEnforcer
     * .resetCounters()} runs on stream creation, so at most a strike or two can accrue against
     * {@code MAX_PING_STRIKES = 2} and they are then reset. Turning it on, by contrast, is precisely
     * what would make an unauthenticated peer unevictable: the server asks for no client certificate
     * and binds every interface (see {@link #LOCAL_HOST}), so a process that completes the TLS
     * handshake and never opens a stream would be free to hold the socket forever.
     *
     * <p>Extracted rather than inlined into the constructor so a spec can assert the four-way choice
     * -- three calls made, one deliberately not made -- against a mock builder. There is no way to
     * read these back off a built {@link Server}, so without the seam a dropped call in a later edit
     * would be invisible to the suite.
     */
    @PackageScope
    static void applyKeepAlive(ServerBuilder<?> builder) {
        // written as separate statements, not a chain: the receiver is wildcard-typed, and Groovy's
        // static checker handles a chain off `ServerBuilder<?>` inconsistently
        builder.keepAliveTime(KEEPALIVE_TIME_SECONDS, TimeUnit.SECONDS)
        builder.keepAliveTimeout(KEEPALIVE_TIMEOUT_SECONDS, TimeUnit.SECONDS)
        builder.permitKeepAliveTime(PERMIT_KEEPALIVE_TIME_SECONDS, TimeUnit.SECONDS)
    }

    private static AgentRpcConfig configFor(Session session) {
        final nav = session?.config?.navigate('agent.rpc')
        final Map opts = nav instanceof Map ? (Map)nav : Collections.emptyMap()
        return new AgentRpcConfig(opts)
    }

    static synchronized AgentRpcBroker get() {
        if( singleton == null ) {
            final session = Global.session as Session
            singleton = new AgentRpcBroker(configFor(session), session)
        }
        return singleton
    }

    @PackageScope
    void close() {
        // Idempotent: `session.onShutdown` and an explicit close (the test harness rebuilds the
        // singleton, and an aborted run can reach both paths) would otherwise report the retention
        // aggregate twice for the same broker.
        if( !closed.compareAndSet(false, true) )
            return
        // Phased shutdown: stop accepting new work and let in-flight dispatch drain
        // before the server is torn down, so a finishing dispatch does not throw in
        // send(). Fall back to a forced shutdown if graceful draining does not complete.
        dispatchPool.shutdown()
        server.shutdown()
        try {
            if( !dispatchPool.awaitTermination(5, TimeUnit.SECONDS) )
                dispatchPool.shutdownNow()
            if( !server.awaitTermination(5, TimeUnit.SECONDS) )
                server.shutdownNow()
        }
        catch( InterruptedException e ) {
            dispatchPool.shutdownNow()
            server.shutdownNow()
            Thread.currentThread().interrupt()
        }
        // before the clear, so the report can count what was still pending
        expiryPool.shutdownNow()
        reportRetention()
        invocations.clear()
        outcomes.clear()
        synchronized(AgentRpcBroker) {
            if( singleton != null && singleton.is(this) )
                singleton = null
        }
    }

    /**
     * Reports, once per run, the capabilities that were minted and never used.
     *
     * <p>This is the only place the shape it exists for is visible. {@code register()} runs while the
     * task SCRIPT is generated -- {@code TaskProcessor} calls {@code task.resolve(taskBody)} BEFORE
     * it consults {@code storeDir} and before the resume cache -- so on {@code nextflow run -resume}
     * every agent task that turns out to be a cache HIT still mints a capability that nothing ever
     * connects with and nothing ever cancels. Counting FIRED deadlines cannot surface that: such a
     * run is over in seconds and the default budget is an hour, so zero deadlines fire and the
     * operator gets silence. What is still PENDING at shutdown is the number that names it.
     *
     * <p>Note what a pending capability holds. The {@link Invocation} keeps the whole
     * {@code AgentRunnerRequest}, including {@code dispatch} -- a {@code ConvertedClosure} over the
     * task-body closure whose {@code delegate} is {@code TaskRun.context} -- so it transitively pins
     * every resolved INPUT VALUE of its task. A wide agent fan-out resumed from cache holds all of
     * it, and widening the pre-connect budget from 180s to an hour widened that window twentyfold.
     */
    private void reportRetention() {
        final pending = invocations.size()
        final lapsed = lapsedCount.get()
        if( !pending && !lapsed )
            return
        final List<String> parts = []
        if( pending )
            parts << "${capabilityCount(pending)} registered but never consumed".toString()
        if( lapsed )
            parts << "${capabilityCount(lapsed)} that lapsed on the ${capabilitySeconds}s `agent.rpc.capabilityTimeout`".toString()
        log.warn "Agent RPC broker shut down with ${parts.join(' and ')} - a capability is minted while the task script is generated, which happens before `storeDir` and the resume cache are consulted, so an agent task that turns out to be a cache hit leaves one behind holding its request, hence that task's resolved inputs, until the run ends or the budget elapses"
    }

    private static String capabilityCount(int count) {
        return count == 1 ? '1 agent capability' : "${count} agent capabilities".toString()
    }

    AgentRpcRegistration register(AgentRunnerRequest request, boolean remote) {
        // the request's own address first: it is the one the guard resolved FOR THIS AGENT, with
        // context the fallback cannot see
        final advertised = remote
                ? (request.brokerHost ?: fallbackHost)
                : AgentRpcHost.of(LOCAL_HOST, 'in-process runner')
        final host = advertised?.host
        if( !host )
            // the ladder's own message names the row that refused, what was tried and what to set;
            // replacing it with a generic line here would discard exactly that
            throw new IllegalStateException(advertised?.error
                    ?: "No address is configured for the driver's agent RPC broker - set `agent.rpc.remoteHost` to a host the agent task can reach the driver on")
        final id = UUID.randomUUID().toString()
        final bytes = new byte[32]
        random.nextBytes(bytes)
        final token = Base64.getUrlEncoder().withoutPadding().encodeToString(bytes)
        // The BROKERED names only. A runner-native tool (`request.nativeToolNames`, the pi SDK
        // builtins the harness adds to the session allowlist) is executed inside the runner
        // container and has no dispatcher on this side, so admitting one here would let the model
        // relocate, say, `bash` from the container into the driver JVM by simply calling it over
        // the RPC stream. brokeredToolNames() is what refuses to build such an allowlist.
        final Set<String> allowedTools = request.brokeredToolNames()
        final invocation = new Invocation(id: id, token: token, request: request, allowedTools: allowedTools)
        invocations.put(id, invocation)
        // A capability's lifetime is the PRE-CONNECT wait, and nothing else. register() runs while
        // the task SCRIPT is generated -- TaskProcessor resolves the body long before it submits the
        // job -- so this clock has to absorb the executor's queueing latency. Deriving it from the
        // per-request LLM timeout (max(requestTimeout,30)+60, i.e. 180s by default) made any queueing
        // executor fail with a security-shaped `Invalid agent RPC invocation identity or token`, so
        // it is now its own generous, configurable budget. Once the connect frame is accepted the
        // capability is consumed and this task is cancelled; the live stream is not on a timer.
        //
        // The deadline is also the ONLY thing that ever releases a capability nobody connects with,
        // and register() runs on paths where no job is ever submitted -- a resume cache hit, a
        // `storeDir` hit -- so those hold their request for the full budget. Three narrower bounds
        // were considered and each is worse than the budget it would replace:
        //
        //  - Revoking on "this task will never run". There is nothing to revoke BY: an
        //    AgentRunnerRequest carries no task identity of any kind (its workDir is a literal '.'),
        //    so a cached TaskRun cannot be mapped back to an invocation id. Building that mapping
        //    means having this plugin observe TaskRun lifecycle across a module boundary -- through
        //    a TraceObserverFactory this plugin does not have -- which is the coupling this design
        //    rejects, and it is only ever an approximation of the terminal state.
        //  - Capping the map. Evicting the OLDEST pending capability is actively harmful: under a
        //    FIFO scheduler the oldest queued tasks are the ones about to start, so the cap would
        //    kill exactly the tasks it must not. Blocking register() at the cap instead deadlocks on
        //    the very shape this is about -- on an all-cached resume nothing ever connects, so the
        //    cap never drains. There is no safe overflow policy.
        //  - Holding the request weakly until connect. The start frame is assembled FROM the request,
        //    and nothing can reconstitute it, so a cleared reference is an unrecoverable failure at
        //    exactly the moment the task finally starts.
        //
        // So the retention stands, bounded by min(capabilityTimeout, run lifetime) -- close() clears
        // the map -- and reportRetention() makes it countable instead of silent.
        invocation.expiry = expiryPool.schedule({ lapse(id, invocation) } as Runnable, capabilitySeconds, TimeUnit.SECONDS)
        // The consumer can only cancel a deadline it can see, and it may have connected -- or the
        // deadline may already have fired, at a capabilityTimeout of zero -- while this one was
        // being armed. Absence from the map is exactly "no longer pending", so re-check and cancel
        // here rather than leave an orphan holding the request for the full budget.
        if( !invocations.containsKey(id) )
            invocation.expiry.cancel(false)
        // the SOURCE is not decoration: an inferred address is indistinguishable from a configured
        // one once it is on the wire, and the expensive failure of the whole ladder is a
        // plausible-but-unroutable one, so name the row that answered -- and shout the rows that
        // resolved something they are not certain about (a multi-homed driver, a containerized
        // driver whose task may land on another docker network)
        final endpoint = endpoint(host, server.port)
        log.debug "Registering agent RPC invocation ${id} remote=${remote} advertising ${endpoint} (${advertised.source})"
        if( advertised.warnings && hostWarned.add(host) )
            for( final warning : advertised.warnings )
                log.warn "Agent RPC broker is advertising ${endpoint} - ${warning}"
        return new AgentRpcRegistration(id, token, endpoint, fingerprint, insecureTransport)
    }

    /**
     * {@code host:port}, with an IPv6 literal BRACKETED. The endpoint travels on argv as
     * {@code --endpoint} and is parsed by Go's {@code grpc.NewClient}, and neither reads a bare
     * {@code 2001:db8::1:8080} as an address and a port. A resolved host can be IPv6 whenever the
     * driver's only routable interface address is (an IPv6-only fabric).
     */
    protected static String endpoint(String host, int port) {
        return host.contains(':') ? "[${host}]:${port}".toString() : "${host}:${port}".toString()
    }

    /**
     * Releases a capability nobody ever connected with, and says so.
     *
     * <p>The previous body was {@code { invocations.remove(id, invocation) }} with no logging at all,
     * which is why widening the budget cut the FREQUENCY of the misleading rejection twentyfold and
     * left its diagnosability at zero: an operator whose queue outran the budget still got
     * {@code Invalid agent RPC invocation identity or token} and not one line of evidence that a
     * deadline had fired. Logging only when the {@code remove} actually SUCCEEDS keeps this quiet for
     * a capability that was consumed a microsecond earlier and lost the race.
     */
    private void lapse(String id, Invocation invocation) {
        if( !invocations.remove(id, invocation) )
            return
        outcomes.put(id, Outcome.EXPIRED)
        final message = "Agent RPC capability ${id} expired after ${capabilitySeconds}s without a connection - the agent task never dialled back within `agent.rpc.capabilityTimeout`"
        // A wide fan-out lapses en masse (a tightened budget, or a run that outlives its queue), and
        // one line per task would bury the run's real errors. The first few carry the diagnosis; the
        // aggregate in reportRetention() carries the true total whatever this cap swallows.
        if( lapsedCount.incrementAndGet() <= LAPSE_WARN_LIMIT )
            log.warn(message)
        else
            log.debug(message)
    }

    /**
     * The PEM of the certificate this broker serves, or {@code null} when TLS is disabled. Exposed
     * for a test client that has to trust an identity no CA vouches for; a certificate is public
     * material, and the private key is never exposed.
     */
    @PackageScope
    String getCertificatePem() { certificatePem }

    private BindableService service() {
        return new BindableService() {
            @Override ServerServiceDefinition bindService() {
                return ServerServiceDefinition.builder(SERVICE_NAME)
                    .addMethod(CONNECT_METHOD, ServerCalls.asyncBidiStreamingCall({ StreamObserver<String> responses ->
                        connect(responses)
                    } as ServerCalls.BidiStreamingMethod<String,String>))
                    .build()
            }
        }
    }

    private StreamObserver<String> connect(StreamObserver<String> responses) {
        final ResponseSink sink = new ResponseSink(responses)
        return new StreamObserver<String>() {
            Invocation invocation
            /** Set once the call has been failed with a status, so frames already on the wire when
             * that happened cannot re-enter the handshake behind it. */
            boolean rejected

            @Override void onNext(String frame) {
                if( rejected )
                    return
                final Object parsed = new JsonSlurper().parseText(frame)
                if( !(parsed instanceof Map) )
                    throw new IllegalArgumentException('Agent RPC frame must be a JSON object')
                final Map msg = (Map)parsed
                if( invocation == null ) {
                    if( msg.type != 'connect' )
                        throw new IllegalArgumentException('First agent RPC frame must be connect')
                    final String id = msg.invocationId?.toString()
                    if( !id ) {
                        // A ConcurrentHashMap forbids a null key, so a connect frame with no identity
                        // used to throw out of onNext and close the call as `UNKNOWN: Application
                        // error processing RPC` -- the opaque answer this handshake exists to stop
                        // giving, reachable before authentication by anyone who can reach the port.
                        rejected = true
                        reject(sink, id, null)
                        return
                    }
                    final candidate = invocations.get(id)
                    if( candidate == null || candidate.token != msg.token?.toString() ) {
                        // A wrong token against a LIVE capability discloses nothing beyond the
                        // generic answer, so only the "no longer pending" case consults the record.
                        rejected = true
                        reject(sink, id, candidate == null ? outcomes.get(id) : null)
                        return
                    }
                    // Consume the capability atomically. A leaked command-line token cannot
                    // reconnect or race a second stream after the legitimate proxy connects. Losing
                    // this CAS means the expiry runnable -- or a second connect -- got there in the
                    // nanoseconds since the get() above, and whichever won recorded WHAT it did, so
                    // report that instead of guessing (this used to report a lapsed capability as a
                    // replay).
                    if( !invocations.remove(candidate.id, candidate) ) {
                        rejected = true
                        reject(sink, id, outcomes.get(id))
                        return
                    }
                    outcomes.put(id, Outcome.CONSUMED)
                    invocation = candidate
                    candidate.expiry?.cancel(false)
                    sink.send(startFrame(candidate))
                    return
                }
                if( msg.invocationId != invocation.id ) {
                    log.warn "Rejected agent RPC frame with mismatched invocation identity (expected=${invocation.id})"
                    rejected = true
                    // same reason as the connect rejections: a thrown exception reaches the task as a
                    // bare UNKNOWN, so the status has to be sent explicitly
                    sink.fail(Status.PERMISSION_DENIED.withDescription('Mismatched agent RPC invocation identity'))
                    return
                }
                switch( msg.type ) {
                    case 'trace':
                        trace(invocation.request, msg)
                        break
                    case 'tool_call':
                        dispatchTool(invocation, msg, sink)
                        break
                    case 'complete':
                        break
                    case 'error':
                        final detail = AgentSecretMasker.redact((msg.message ?: msg.code)?.toString(), invocation.request.apiKey)
                        log.warn "Agent `${invocation.request.agentName}` failed remotely: ${detail}"
                        break
                    default:
                        throw new IllegalArgumentException("Unknown agent RPC frame type: ${msg.type}")
                }
            }

            @Override void onError(Throwable error) {
                sink.markClosed()
                if( invocation != null )
                    log.debug("Agent RPC stream failed for ${invocation.id}", error)
            }

            @Override void onCompleted() {
                sink.complete()
            }
        }
    }

    /**
     * The {@code start} frame: the portable {@link AgentProtocolSpec}, and BESIDE it the provider
     * credential core resolved for this request. Beside and never inside -- the spec is the payload
     * a transport may relay, persist or log, and it stays credential-free by construction.
     *
     * <p>The frame carries {@code request.apiKey} and NEVER {@code request.credential()}. The
     * placeholder that method substitutes when an endpoint is declared but nothing resolved is a
     * driver-side artifact of a client that refuses to issue a request without a key; a runtime key
     * OWNS its provider in pi ({@code setRuntimeApiKey} is consulted before pi's own store and
     * before the ambient environment), so sending {@code nxf-no-credential} here would shadow the
     * very credential an {@code env}/{@code secret} channel had delivered to the container and turn
     * a working run into a 401.
     *
     * <p>Sending the real one is what naming the provider makes safe. The broker used to withhold it
     * because a runtime key owns its provider and the driver could not establish that the key
     * belonged to THIS provider; {@link nextflow.agent.AgentConfig#apiKeyFor} now answers that -- it
     * resolves in the model's own provider namespace and refuses to hand over an ambient provider
     * variable when the endpoint belongs to somebody else -- so what arrives here is already scoped
     * to the provider the task will call.
     *
     * <p>Gated on transport security. Under {@code agent.rpc.tls = false} the frame is cleartext and
     * the driver is not authenticated to the task, so the credential is withheld and the run is told
     * how to deliver it out of band: the escape hatch keeps working, it just carries no secret.
     *
     * <p>NOTE the exposure this accepts even WITH TLS. The capability token is on the task's argv,
     * hence in {@code .command.sh} in the work directory, so from here on read access to a work
     * directory is read access to the provider credential -- for as long as the capability is
     * pending ({@code agent.rpc.capabilityTimeout}, an hour by default, and a resume cache hit mints
     * capabilities nothing ever consumes). Single use and a pinned certificate bound a stolen
     * token to ONE connection -- not to none.
     */
    private Map startFrame(Invocation invocation) {
        final Map frame = [
            type: 'start',
            protocolVersion: 2,
            invocationId: invocation.id,
            spec: AgentProtocolSpec.fromRequest(invocation.request) ] as Map
        final credential = invocation.request.apiKey
        if( invocation.request.credentialWithheld )
            warnWithheldCredential(invocation.request)
        if( !credential )
            return frame
        if( insecureTransport ) {
            warnInsecureCredential()
            return frame
        }
        frame.put('apiKey', credential)
        return frame
    }

    /**
     * A provider credential resolved on the driver and the endpoint gate refused to send it
     * ({@code AgentConfig.credentialWithheldFor}). WARN, never throw.
     *
     * <p>This asymmetry with the langchain4j runner -- which raises the same condition as a fatal
     * error in {@code ChatModelFactory.createModel} -- is deliberate and must stay. langchain4j is
     * the whole credential chain: core resolving nothing means nothing exists. pi is not: it reads
     * its own auth store and the provider variables present in the CONTAINER, which the driver
     * cannot see and which are exactly how a Kubernetes Secret or an {@code env}/{@code secret}
     * channel delivers a key. Aborting here would break a deployment that is working.
     *
     * <p>Said once per model id: the cause is that model's configuration, so a fan-out of tasks
     * would otherwise restate one fact hundreds of times.
     */
    private void warnWithheldCredential(AgentRunnerRequest request) {
        if( !withheldCredentialWarned.add(request.model ?: '') )
            return
        log.warn "The LLM provider credential resolved for agent model `${request.model}` was withheld from ${request.baseUrl ? "the endpoint ${request.baseUrl}" : 'the provider default endpoint'} because that endpoint is not one the `${request.apiProvider}` namespace owns - it is NOT sent to the agent task; the runner falls back to its own credential store and to the provider variables in the container, or set `agent.apiKey` to the credential this endpoint accepts, or `agent.apiProvider` to name the namespace it belongs to"
    }

    /**
     * Says why the credential was withheld, ONCE per broker: the cause is the session's
     * {@code agent.rpc.tls}, not this invocation, so a line per agent task would restate one fact a
     * fan-out's worth of times over the run's real errors.
     */
    private void warnInsecureCredential() {
        if( !insecureCredentialWarned.compareAndSet(false, true) )
            return
        log.warn "Agent RPC transport security is disabled (agent.rpc.tls=false) - the LLM provider credential resolved by the driver is NOT sent to the agent task, because the start frame is cleartext and the driver is not authenticated to it; re-enable `agent.rpc.tls`, or deliver the credential to the container out of band with `agent.containerOptions = '-e OPENAI_API_KEY'` (Docker/Podman), the `env` config scope, or the `secret` directive"
    }

    /**
     * Refuses a {@code connect} frame, telling the task WHICH of the three ways it can fail happened.
     *
     * <p>{@code outcome} is what the pending map no longer holds: {@code null} means the id was never
     * issued (or its token did not match a live capability), which is the only case that warrants the
     * indiscriminate answer. The other two used to produce that same answer, and that is the whole
     * defect -- a queue longer than the budget, and a task retried after it had already connected,
     * are both operational events being reported as if they were forged credentials.
     *
     * <p>The record is written just AFTER the winning {@code remove}, never before, so that a lapse
     * which loses the CAS cannot stamp {@code EXPIRED} over a capability that was in fact consumed.
     * The cost is a nanosecond-wide window in which a connect sees the entry gone but no outcome yet
     * and falls back to the generic answer -- a degraded message, never a wrong one.
     */
    private void reject(ResponseSink sink, String id, Outcome outcome) {
        String detail
        if( outcome == Outcome.EXPIRED ) {
            detail = "Agent RPC invocation capability expired after ${capabilitySeconds}s while the task waited to start - raise `agent.rpc.capabilityTimeout` above the executor's queueing delay".toString()
            log.warn "Rejected agent RPC connection with a capability that expired after ${capabilitySeconds}s (invocationId=${id})"
        }
        else if( outcome == Outcome.CONSUMED ) {
            // Named cause first, replay second, and on purpose. TaskProcessor re-submits a
            // ProcessRetryableException / CloudSpotTerminationException copy via task.makeCopy()
            // WITHOUT re-resolving the task body -- unlike the ordinary retry path, which does -- so a
            // node termination or a spot reclaim re-runs the SAME `--invocation`/`--token`. If the
            // first attempt had got as far as connecting, that legitimate retry lands here, and
            // labelling it a replay attack is precisely the security-shaped mislabel this PR removes.
            detail = 'Agent RPC invocation capability was already consumed - the task was retried after an earlier attempt had already connected, or the capability was replayed'
            log.warn "Rejected agent RPC connection presenting an already-consumed capability (invocationId=${id}) - usually a task retried after an earlier attempt had connected; a replayed capability is the other reading"
        }
        else {
            detail = 'Invalid agent RPC invocation identity or token'
            // Never log the token value; the claimed invocation id is enough to audit.
            log.warn "Rejected agent RPC connection with invalid identity or capability token (invocationId=${id})"
        }
        sink.fail(Status.UNAUTHENTICATED.withDescription(detail))
    }

    private void dispatchTool(Invocation invocation, Map msg, ResponseSink sink) {
        dispatchPool.submit {
            final String callId = msg.callId?.toString()
            final String name = msg.name?.toString()
            try {
                if( !callId || !name )
                    throw new IllegalArgumentException('Invalid tool_call without callId/name')
                if( !invocation.allowedTools.contains(name) )
                    throw new SecurityException("Tool `${name}` is not authorized for this invocation")
                if( !invocation.callIds.add(callId) )
                    throw new IllegalArgumentException("Duplicate tool call identity `${callId}`")
                // Enforce the call ceiling atomically: two dispatch threads could otherwise
                // both observe a within-limit size between add() and the size() check.
                if( invocation.callCount.incrementAndGet() > Math.max(invocation.request.maxIterations, 1) ) {
                    invocation.callCount.decrementAndGet()
                    invocation.callIds.remove(callId)
                    throw new IllegalStateException('Agent exceeded the authorized tool-call limit')
                }
                final dispatcher = invocation.request.dispatch
                if( dispatcher == null )
                    throw new IllegalStateException("No dispatcher is available for tool `${name}`")
                final args = JsonOutput.toJson(msg.arguments ?: Collections.emptyMap())
                final result = dispatcher.call(name, args)
                trySend(sink, invocation, callId, [
                    type: 'tool_result',
                    invocationId: invocation.id,
                    callId: callId,
                    result: result,
                    isError: false ])
            }
            catch( Throwable error ) {
                log.warn("Agent RPC tool dispatch failed for invocation ${invocation.id} tool=${name} callId=${callId}", error)
                trySend(sink, invocation, callId, [
                    type: 'tool_result',
                    invocationId: invocation.id,
                    callId: callId,
                    result: error.message ?: error.toString(),
                    isError: true ])
            }
        }
    }

    /**
     * Sends a frame through the response sink, containing any failure so a dead
     * stream cannot abort the dispatch task or escape into the never-read Future.
     */
    private static void trySend(ResponseSink sink, Invocation invocation, String callId, Map message) {
        try {
            sink.send(message)
        }
        catch( Throwable error ) {
            log.warn("Failed to send agent RPC frame for invocation ${invocation.id} callId=${callId}", error)
        }
    }

    private static void trace(AgentRunnerRequest request, Map msg) {
        if( !request.trace )
            return
        final event = msg.event ?: 'event'
        // Trace frames can restate the prompt or model thinking; keep them at DEBUG
        // and redact any embedded credential before logging.
        final detail = AgentSecretMasker.redact((msg.name ?: msg.text ?: '').toString(), request.apiKey)
        log.debug "[agent:${request.agentName ?: 'agent'}] ${event}${detail ? ' ' + detail : ''}"
    }
}
