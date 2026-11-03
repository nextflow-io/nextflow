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

package nextflow.agent.rpc

import java.lang.ref.WeakReference
import java.util.regex.Pattern
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.container.SmolVmConfig
import nextflow.executor.AbstractGridExecutor
import nextflow.executor.Executor
import nextflow.agent.AgentConfig

/**
 * The single owner of the ladder that decides which address a containerized agent task uses to
 * reach the driver's RPC broker.
 *
 * <p>Every caller -- the pre-ignition guard
 * {@link nextflow.agent.AgentLaunchConditions#requireBrokerHost}, the broker that advertises the
 * endpoint, and {@link AgentRpcConfig#resolveRemoteHost} -- goes through here, so there is one
 * implementation to keep in agreement rather than three. The rows are evaluated in STRICT order,
 * first hit wins:
 *
 * <pre>
 *  R1  agent.rpc.remoteHost set                            -> that value
 *  R2  NXF_AGENT_RPC_REMOTE_HOST set                       -> that value
 *  R3  local + docker/podman + local daemon + containerized driver (not rootless) -> outbound address
 *  R4  local + the container shares the host network namespace, and that namespace
 *      belongs to THIS kernel (not to a Docker Desktop / podman machine VM)       -> 127.0.0.1
 *  R5  local + docker/podman + local daemon                -> the engine host alias
 *  R6  local + a network-isolated engine with no alias     -> outbound address
 *  R7  the executor is an AbstractGridExecutor             -> outbound address
 *  R8  k8s, with the client config resolved FROM the cluster                      -> outbound address
 *  R9  a cloud batch executor, with a positive cloud-membership probe             -> outbound address
 *  R10 otherwise                                           -> error row E7
 * </pre>
 *
 * <p>The "outbound address" is the local address the kernel picks for the default route -- a
 * routing-table lookup, no packet sent. It was measured reachable from every container engine that
 * could be tested on the probe host, including the two VM-isolated ones, which is why a single
 * address is advertised and never a candidate list: the capability token may only be offered to a
 * driver whose TLS fingerprint has already verified, so trying N candidates would offer it to N
 * hosts.
 *
 * <p>The cases that CANNOT work are error rows E1..E7 rather than a guess. Each carries its own
 * message naming what was tried and what to set, and each is raised BEFORE ignition: a capability
 * is only ever released by its one-hour {@code agent.rpc.capabilityTimeout}, so letting a doomed
 * configuration submit trades a loud failure for an hour-long silent one.
 *
 * <p>Everything the ladder observes about the host -- the routing table, whether the driver is
 * itself containerized, whether it runs rootless, whether it is a cloud instance -- comes through
 * {@link Probes}, injected so every row is unit-testable with no network, no container and no
 * cluster, and memoized so the probes run at most ONCE per session however many agent definitions
 * the script declares.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AgentRpcHostResolver {

    /** The address a container in the driver's own network namespace reaches the driver on. */
    static final String LOOPBACK = '127.0.0.1'

    /** Engines with a daemon of their own, hence an alias for the daemon host and a remote mode. */
    private static final List<String> DAEMON_ENGINES = ['docker', 'podman']

    /**
     * Engines that create NO network namespace by default -- the container simply shares the
     * driver's -- but which accept a run option that makes them create one. Inspected in the
     * OPPOSITE direction from docker/podman: the option is what disqualifies the row (E6).
     */
    private static final String SMOLVM = 'smolvm'

    private static final String APPLE_CONTAINER = 'apple-container'

    /** Engines that create no network namespace at all, with no switch to make them. */
    private static final List<String> HOST_NAMESPACE_ENGINES = ['apptainer', 'singularity']

    /** Executors that run the agent task on an instance of a cloud compute environment. */


    /**
     * The conventional default-bridge range per engine. A containerized driver whose own address
     * falls OUTSIDE it is very likely on a user-defined network (a {@code docker compose} run, or
     * {@code docker run --network <name>}), while Nextflow emits no {@code --network} for the agent
     * task, so the task lands on the default bridge -- and docker's own DOCKER-ISOLATION rules DROP
     * traffic between bridges. That is engine behaviour, not a site firewall, so it earns a warning.
     */
    private static final Map<String,String> DEFAULT_BRIDGE_PREFIX = [docker: '172.17.', podman: '10.88.']


    private final Map sessionConfig
    private final Probes probes

    // -- memoized probe results; each is computed at most once for the life of the resolver, which
    //    is at most once per session (see #of). `outboundResolved` distinguishes "not looked up yet"
    //    from "looked up, and the host has no default route".
    private boolean outboundResolved
    private String outboundValue
    private List<String> interfacesValue
    private Boolean containerizedValue
    private Boolean rootlessValue
    private Boolean cloudValue

    AgentRpcHostResolver(Map sessionConfig, Probes probes) {
        this.sessionConfig = sessionConfig
        this.probes = probes
    }

    // ---------------------------------------------------------------------------------------
    // per-session instance
    // ---------------------------------------------------------------------------------------

    private static final Object LOCK = new Object()
    /**
     * WEAK, deliberately: this is a static field on a class that lives for the JVM, so a strong
     * reference would pin the last {@link Session} -- and transitively its whole config map -- past
     * the end of its run. An embedded host that executes several pipelines in one JVM would retain
     * every session but the current one. Identity is all the cache key needs, and a cleared
     * reference simply misses.
     */
    private static WeakReference<Session> cachedSession
    private static AgentRpcHostResolver cachedResolver

    /**
     * The resolver for this session, created once. Memoization has to outlive a single call because
     * {@code AgentLaunchConditions.requireCanonicalLaunch} runs once per agent DEFINITION while
     * every probe answers a question about the HOST, so a script with ten agents would otherwise
     * open ten sockets and read {@code /proc} ten times for ten identical answers.
     *
     * <p>The {@link AgentRpcConfig} is deliberately NOT part of the cached state: it is rebuilt per
     * agent definition ({@code AgentDef.agentConfig()}) and a selector can give two agents different
     * {@code agent.rpc.remoteHost} values, so it is passed per call while only the host facts -- the
     * ones this cache exists for -- are shared. With no session there is nothing to key on, so each
     * call gets its own resolver rather than inheriting another caller's memo.
     */
    static AgentRpcHostResolver of(Session session) {
        if( session == null )
            return new AgentRpcHostResolver(null, new SystemProbes())
        synchronized (LOCK) {
            if( cachedResolver == null || cachedSession?.get() !== session ) {
                cachedSession = new WeakReference<Session>(session)
                cachedResolver = new AgentRpcHostResolver(session.getConfig(), new SystemProbes())
            }
            return cachedResolver
        }
    }

    /**
     * Resolve the address for one agent definition. The {@code executor} INSTANCE and its NAME are
     * both required and neither is redundant: R7 keys off {@code instanceof AbstractGridExecutor},
     * because the executor registry carries names a hand-written list would miss, while R8/R9 key
     * off the {@code @ServiceName} string.
     */
    static AgentRpcHost resolve(Executor executor, String executorName, ContainerConfig containerConfig,
            Object containerOptions, AgentRpcConfig rpc, Session session) {
        return of(session).resolve(executor, executorName, containerConfig, containerOptions, rpc)
    }

    /** Discard the per-session instance; for tests, which run many sessions in one JVM. */
    static void reset() {
        synchronized (LOCK) {
            cachedSession = null
            cachedResolver = null
        }
    }

    /**
     * Seed the per-session instance with alternative {@link Probes}, so a test can drive the ladder
     * through the callers that reach it statically -- the pre-ignition guard and the broker -- and
     * still assert which RUNG answered rather than which machine the suite happens to run on. Without
     * it a suite running in a container would see {@code /.dockerenv} and take R3 for every alias row.
     * Paired with {@link #reset}.
     */
    static AgentRpcHostResolver install(Session session, Probes probes) {
        synchronized (LOCK) {
            cachedSession = session != null ? new WeakReference<Session>(session) : null
            cachedResolver = new AgentRpcHostResolver(session?.getConfig(), probes)
            return cachedResolver
        }
    }

    // ---------------------------------------------------------------------------------------
    // the ladder
    // ---------------------------------------------------------------------------------------

    /**
     * Whether the container engine that launches the agent task runs on the DRIVER host, which is a
     * NECESSARY but not sufficient condition for R3..R6: it proves the EXECUTOR is local, not that
     * the container daemon is (see E3).
     */
    static boolean isDriverHostExecutor(String executor) {
        return executor == null || executor == AgentConfig.DEFAULT_EXECUTOR
    }

    AgentRpcHost resolve(Executor executor, String executorName, ContainerConfig containerConfig,
            Object containerOptions, AgentRpcConfig rpc) {
        return resolve(executor, executorName, containerConfig?.getEngine(), containerConfig, containerOptions, rpc)
    }

    /**
     * The engine NAME is taken separately from the {@link ContainerConfig} so a caller that has only
     * the name -- {@link AgentRpcConfig#resolveRemoteHost}, which is the shape the runner SPI has
     * always had -- goes through the same ladder rather than around it. Such a caller gets no engine
     * {@code runOptions} and no {@code smolvm.network}, so it can only be answered by the rows that
     * do not need them; that is a property of what it knows, not a second ladder.
     *
     * <p>{@code synchronized} because the memo fields below are written here: the pre-ignition guard
     * calls this from the SCRIPT thread while the broker's fallback ({@link
     * AgentRpcConfig#resolveBrokerHost}) can call it from a task thread, and there is no other
     * happens-before edge between the two. The probes are idempotent, so the lock only ever costs a
     * repeated probe it also makes unnecessary.
     */
    synchronized AgentRpcHost resolve(Executor executor, String executorName, String engine, ContainerConfig containerConfig,
            Object containerOptions, AgentRpcConfig rpc) {
        return ladder(executor, executorName, engine, containerConfig, containerOptions, rpc)
    }

    private AgentRpcHost ladder(Executor executor, String executorName, String engine,
            ContainerConfig containerConfig, Object containerOptions, AgentRpcConfig rpc) {
        final List<String> tried = []

        // -- R1/R2: an explicit value, from the config or from the environment. AgentRpcConfig has
        //    already merged the two rungs, and remembers which one answered so the log line can say.
        if( rpc?.getRemoteHost() )
            return AgentRpcHost.of(rpc.getRemoteHost(), rpc.getRemoteHostSource())
        tried << '`agent.rpc.remoteHost` (not set)'
        tried << "`${AgentRpcConfig.REMOTE_HOST_ENV}` (not set)".toString()

        // -- The engine decides the answer ONLY when it runs on the driver's own machine. Anywhere
        //    else, where the task runs is the executor's business and the engine knows nothing
        //    about how the task reaches back.
        if( isDriverHostExecutor(executorName) ) {
            // -- R3: docker and podman synthesize a name for the host their containers run on.
            if( engine in DAEMON_ENGINES )
                return AgentRpcHost.of(AgentRpcConfig.hostAliasFor(engine), "${engine} host alias".toString())

            // -- R4: these engines create no network namespace, so the driver's loopback IS the
            //    container's loopback. A property of the engine, not a guess about the deployment.
            if( engine in HOST_NAMESPACE_ENGINES )
                return AgentRpcHost.of(LOOPBACK, 'host network namespace')

            // -- R5: a microVM has its own kernel, so loopback is the GUEST's. With networking
            //    disabled there is no address at all and no rung below can help.
            if( engine == SMOLVM ) {
                if( !smolvmNetwork(containerConfig) )
                    return smolvmNetworkError()
                return outbound('inferred from default route', tried)
            }

            // -- R6: apple-container is VM-isolated and synthesizes no host name of its own. The
            //    outbound address was measured reachable from it.
            if( engine == APPLE_CONTAINER )
                return outbound('inferred from default route', tried)

            // No engine at all means nothing is containerized here, so there is no task to reach
            // back and no address to infer. An UNKNOWN engine is different: it containerizes, it
            // simply has no rung of its own, so it takes the routable address below like any
            // remote deployment.
            if( !engine ) {
                tried << 'the container engine (none is enabled, so no agent task is containerized)'
                return unresolvedError(tried)
            }
            tried << "the container engine (`${engine}` names no address for the driver host)".toString()
        }

        // -- R7: everything else -- Kubernetes, a grid/HPC cluster, a cloud batch service, or a
        //    driver-host engine with no rung of its own. They share one shape: the task runs
        //    somewhere that must route back to the driver. Deliberately NOT split per executor:
        //    the answer is identical for all of them, so a per-provider rung would be a branch
        //    with no different behaviour behind it, and a provider probe we cannot test.
        return outbound('inferred from default route', tried)
    }

    // ---------------------------------------------------------------------------------------
    // the outbound address
    // ---------------------------------------------------------------------------------------

    /**
     * The default-route address, checked for usability. A plausible-but-unroutable address is the
     * expensive failure this whole design exists to avoid, so loopback, wildcard, link-local and
     * "no address at all" are error row E7 rather than something to advertise.
     */
    private AgentRpcHost outbound(String source, List<String> tried, List<String> warnings = null) {
        String address = outboundAddress()
        String label = source
        if( !usableAddress(address) ) {
            // A host with NO default route is not necessarily offline: an air-gapped cluster and an
            // IPv6-only fabric both route perfectly well on the interface the compute nodes share,
            // and that address is already enumerated one line down. Only an UNAMBIGUOUS answer is
            // taken -- with several usable addresses there is no evidence which one the agent task
            // can reach, and advertising a guess is the failure this design exists to avoid.
            final candidates = interfaceAddresses().findAll { usableAddress(it) }
            // A non-local deployment may legitimately have to be reached on a PUBLIC address, so a
            // single globally-routable candidate settles an otherwise ambiguous set. Deliberately
            // read off the local interfaces rather than asked of a provider metadata service: no
            // per-cloud probe, nothing to test against, and a NAT-assigned public address is not
            // on the interface anyway -- for which `agent.rpc.remoteHost` is the answer.
            final publics = candidates.findAll { publicAddress(it) }
            if( publics.size() == 1 ) {
                final all = warnings != null ? new ArrayList<String>(warnings) : new ArrayList<String>()
                all << "the driver host has no usable default route, so its public address `${publics[0]}` was advertised - set `agent.rpc.remoteHost` if the agent tasks reach the driver on a different one".toString()
                return AgentRpcHost.of(publics[0], "${source}; the host has no default route, so its only public interface address was used".toString(), all)
            }
            if( candidates.size() != 1 ) {
                tried << "the host's default route (the routing table returned ${address ? '`' + address + '`, which is not routable from another host' : 'no address at all'})".toString()
                tried << (candidates
                        ? "the host's interfaces (${candidates.collect { '`' + it + '`' }.join(', ')} are all routable, so none of them is unambiguously the one to advertise)".toString()
                        : 'the host\'s interfaces (none carries a routable address, so the driver appears to be offline)')
                return unresolvedError(tried)
            }
            address = candidates[0]
            label = "${source}; the host has no default route, so its only routable interface address was used".toString()
        }
        final List<String> all = warnings != null ? new ArrayList<String>(warnings) : new ArrayList<String>()
        // Same ADDRESS FAMILY only. A dual-stack NIC reports its IPv4 and its global IPv6 address
        // separately, so comparing across families would fire this on essentially every ordinary
        // Linux driver and train operators to ignore the one warning §4 makes load-bearing.
        final others = interfaceAddresses().findAll { it != address && sameFamily(it, address) }
        if( others )
            all << "the driver host is multi-homed: the default route selected `${address}`, but ${others.collect { '`' + it + '`' }.join(', ')} also exist - set `agent.rpc.remoteHost` if the agent tasks reach the driver on one of those instead".toString()
        return AgentRpcHost.of(address, label, all)
    }

    /**
     * Whether two address LITERALS belong to the same family. Every value here comes from
     * {@code InetAddress.getHostAddress}, which spells IPv6 with colons and IPv4 without and never
     * emits the mapped {@code ::ffff:} form, so the test needs no parsing and no lookup.
     */
    protected static boolean sameFamily(String a, String b) {
        return a?.contains(':') == b?.contains(':')
    }

    /**
     * The address must be usable FROM ANOTHER HOST. {@code InetAddress} answers all three questions
     * without a lookup, because every value here is already a literal from {@code getHostAddress}.
     */
    protected static boolean usableAddress(String address) {
        if( !address )
            return false
        try {
            final inet = InetAddress.getByName(address)
            return !inet.isLoopbackAddress() && !inet.isAnyLocalAddress() && !inet.isLinkLocalAddress()
        }
        catch( Exception e ) {
            log.debug "Unable to interpret `${address}` as an address for the driver's agent RPC broker - ${e.message}"
            return false
        }
    }

    /**
     * Whether an address is routable from OUTSIDE the driver's own network, i.e. not RFC1918 /
     * unique-local. {@code isSiteLocalAddress} covers 10/8, 172.16/12 and 192.168/16 for IPv4 and
     * {@code fec0::/10} for IPv6; {@code fc00::/7} unique-local is checked separately because
     * {@link InetAddress} has no predicate for it.
     */
    protected static boolean publicAddress(String address) {
        if( !usableAddress(address) )
            return false
        try {
            final inet = InetAddress.getByName(address)
            if( inet.isSiteLocalAddress() )
                return false
            final b = inet.address
            // fc00::/7 - unique local, the IPv6 analogue of RFC1918
            return !(b.length == 16 && (b[0] & 0xFE) == 0xFC)
        }
        catch( Exception e ) {
            return false
        }
    }

    private String outboundAddress() {
        if( !outboundResolved ) {
            outboundValue = probes.outboundAddress()
            outboundResolved = true
        }
        return outboundValue
    }

    private List<String> interfaceAddresses() {
        if( interfacesValue == null )
            interfacesValue = probes.interfaceAddresses() ?: Collections.<String>emptyList()
        return interfacesValue
    }

    // ---------------------------------------------------------------------------------------
    // row conditions
    // ---------------------------------------------------------------------------------------

    /** An IPv4 dotted quad, the only literal shape that needs a pattern to be told from a NAME. */
    private static final Pattern IPV4_LITERAL = Pattern.compile(/^\d{1,3}(?:\.\d{1,3}){3}$/)

    /** {@code --network host}, {@code --network=host}, {@code --net host}, {@code --net=host}. */
    private static final Pattern HOST_NETWORK = Pattern.compile(/(?:^|\s)--net(?:work)?[= ]host(?:\s|$)/)

    /**
     * {@code smolvm.network} defaults to true and {@code SmolVmBuilder} emits {@code --net} only
     * when it is set, so a microVM created with it false has no routes at all.
     */
    protected static boolean smolvmNetwork(ContainerConfig containerConfig) {
        // Read the typed field off the smolvm config rather than reflecting a property name off
        // the generic interface: `ContainerConfig` does not declare `network` (it is engine
        // specific), and anything that is not a SmolVmConfig gets SmolVmBuilder's default, true.
        return !(containerConfig instanceof SmolVmConfig) || ((SmolVmConfig) containerConfig).network
    }

    // ---------------------------------------------------------------------------------------
    // error rows -- one message per row, each naming what was tried and what to set
    // ---------------------------------------------------------------------------------------

    private static AgentRpcHost smolvmNetworkError() {
        return AgentRpcHost.error('E2', 'the `smolvm` microVM is created with no network at all (`smolvm.network = false`), so the agent task cannot dial the driver on any address - set `smolvm.network = true`')
    }

    private static AgentRpcHost unresolvedError(List<String> tried) {
        return AgentRpcHost.error('E7', "no address the agent task could reach the driver on could be determined. Tried, in order: ${tried.join('; ')}. Set `agent.rpc.remoteHost` to a host the agent task can reach the driver on - a routable driver address, an in-cluster service name, or `${LOOPBACK}` when the container shares the driver's network namespace".toString())
    }

    // ---------------------------------------------------------------------------------------
    // the real probes
    // ---------------------------------------------------------------------------------------

}
