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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.exception.ScriptRuntimeException
import nextflow.executor.Executor
import nextflow.agent.rpc.AgentRpcHost
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.rpc.AgentRpcConfig

/**
 * The container and executor conditions a canonical agent launch depends on: whether the resolved
 * configuration will actually run the agent task in a container, which address the task can dial
 * the driver's RPC broker on, and which run options the driver has to add for it.
 *
 * <p>Pure statics that read no agent state -- they change when container engines and executors
 * change, not when agents do -- so they live beside the rest of the agent runtime rather than in
 * {@link nextflow.script.AgentDef}, which merely calls them while lowering an agent to a task.
 *
 * <p>Only three of them are that caller's surface: {@link #requireCanonicalLaunch},
 * {@link #withDockerHostGateway} and {@link #requireTaskContainer}. The rest are the steps those
 * three are composed of and stay {@code @PackageScope} -- visible to this package's tests, which
 * exercise each rung directly, but not to {@code nextflow.script}, which must go through the
 * composed guards so it cannot admit a launch on a subset of the conditions.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AgentLaunchConditions {

    private AgentLaunchConditions() {}

    /** A configured container image excludes both an absent value and the explicit opt-out. */
    @PackageScope
    static boolean hasContainer(Object container) {
        return container != null && container != false
    }

    /**
     * Whether the resolved configuration containerizes the agent task: an image PLUS either an
     * executor that manages containers itself or an enabled container engine. This mirrors
     * {@link nextflow.processor.TaskRun#isContainerEnabled} except that the engine is the one the
     * EXECUTOR asks for ({@link nextflow.executor.Executor#containerConfigEngine}), which is also
     * what the task itself will use -- the engine-agnostic {@code session.getContainerConfig()}
     * disagrees with the task's own view for a container-native executor.
     */
    @PackageScope
    static boolean willContainerize(Object container, boolean containerNative, boolean engineEnabled) {
        return hasContainer(container) && (containerNative || engineEnabled)
    }

    /**
     * Resolve the container facts a canonical agent launch depends on -- whether the executor manages
     * containers itself, and whether the engine that executor asks for is enabled -- and reject the
     * run when they do not add up.
     *
     * <p>The resolved {@link nextflow.executor.Executor} INSTANCE is the oracle, not its name:
     * {@link nextflow.executor.Executor#isContainerNative} can depend on the session (the local
     * executor is container-native under Fusion) and
     * {@link nextflow.executor.Executor#containerConfigEngine} decides WHICH engine block must be
     * enabled. It is the same cached instance {@code createTaskProcessorResolved} reuses. The
     * executor is resolved only once an image is present, so a missing {@code agent.container} is
     * still reported for an executor that cannot be instantiated here.
     *
     * <p>Kept out of {@link nextflow.script.AgentDef#buildAgentTask} so this resolution -- which
     * executor is interrogated and which engine block is read -- is directly unit-testable.
     *
     * @return the container engine that will run the agent task ON THE DRIVER HOST -- or {@code null}
     *      when the container is launched elsewhere (or no image is declared) and the driver's run
     *      options therefore do not apply, see {@link #withDockerHostGateway} -- together with the
     *      resolved broker address, which the caller needs to decide those same run options
     */
    static CanonicalLaunch requireCanonicalLaunch(String agentName, String runner, String executor,
            Object container, Object containerOptions, AgentRpcConfig rpc, Session session,
            boolean runnerAutoSelected = false) {
        boolean containerNative = false
        boolean engineEnabled = false
        String containerEngine = null
        // hoisted out of the block below because the ladder needs BOTH: the executor INSTANCE for
        // its `instanceof AbstractGridExecutor` row, and the engine config for its run options. They
        // are resolved only when an image is present, which is safe only because requireContainerized
        // throws first when it is absent
        Executor agentExecutor = null
        ContainerConfig containerConfig = null
        if( hasContainer(container) ) {
            agentExecutor = session.getExecutorFactory().getExecutorByName(executor, session)
            containerNative = agentExecutor.isContainerNative()
            containerConfig = session.getContainerConfig(agentExecutor.containerConfigEngine())
            engineEnabled = containerConfig != null && containerConfig.isEnabled()
            containerEngine = containerConfig?.getEngine()
        }
        requireContainerized(agentName, runner, executor, container, containerNative, engineEnabled, runnerAutoSelected)
        final brokerHost = requireBrokerHost(agentName, executor, agentExecutor, containerConfig, containerOptions, rpc, session)
        return new CanonicalLaunch(isDriverHostEngine(executor) ? containerEngine : null, brokerHost)
    }

    /**
     * What {@link #requireCanonicalLaunch} establishes: which engine runs the agent task on the
     * driver host, and which address the task will dial the driver's broker on. The two travel
     * together because the run options depend on both -- the docker host-gateway mapping is needed
     * only when the RESOLVED address is the docker alias.
     */
    static class CanonicalLaunch {
        final String containerEngine
        final AgentRpcHost brokerHost

        CanonicalLaunch(String containerEngine, AgentRpcHost brokerHost) {
            this.containerEngine = containerEngine
            this.brokerHost = brokerHost
        }
    }

    /**
     * Whether the container engine that launches the agent task runs on the DRIVER host, which is
     * the only case in which an engine host alias names the driver. True for the local executor
     * (whatever engine it drives, and whether or not Fusion makes it container-native) and false for
     * every other executor: a grid executor's docker daemon runs on a compute node, so its
     * {@code host.docker.internal} names that node's host, not the driver.
     *
     * <p>Delegates to the resolver so the guard and the ladder cannot disagree about what "local"
     * means. It is a NECESSARY but not sufficient condition there: the executor name says nothing
     * about a docker CLI pointed at another machine, which is error row E3.
     */
    @PackageScope
    static boolean isDriverHostEngine(String executor) {
        return AgentRpcHostResolver.isDriverHostExecutor(executor)
    }

    /**
     * Reject a canonical agent whose resolved configuration would NOT run the task in a container:
     * the launch command is built from paths that exist only inside the runner image, so a
     * non-containerized task would fail with `No such file`. Generalises the former
     * offload-only check to every executor -- a grid executor with no container engine is
     * rejected as loudly as a missing image.
     */
    @PackageScope
    static void requireContainerized(String agentName, String runner, String executor,
            Object container, boolean containerNative, boolean engineEnabled, boolean runnerAutoSelected = false) {
        if( willContainerize(container, containerNative, engineEnabled) )
            return
        final hint = runnerAutoSelected ? autoSelectedRunnerHint(runner) : ''
        // reachable only for a runner that has no image of its own ({@link AgentRunner#getDefaultContainer}
        // is null) or for the explicit `agent.container = false` opt-out: a runner that DOES declare
        // an image has already had it defaulted into the config by AgentDef.resolveLaunch
        if( !hasContainer(container) )
            throw new ScriptRuntimeException("Agent `${agentName}` must declare a container - set `agent.container` to the `${runner}` runner image, which carries the agent proxy and harness the task is launched from${hint}")
        // the image may have come from the runner rather than from the user, so this must NOT say
        // the user declared `agent.container` -- for the runner this feature exists for, an image
        // is always present and this is the message the most common first run gets
        throw new ScriptRuntimeException("Agent `${agentName}` has a container image (`agent.container`, or the `${runner}` runner's own image when that is unset) but executor `${executor}` would not run it in a container - enable a container engine (e.g. `docker.enabled = true`) so the `${runner}` runner image is used to run the agent task${hint}")
    }

    /**
     * Names the runner the message just blamed, for the case where the user never chose it:
     * {@link nextflow.agent.AgentRunnerProvider#get} selects the sole installed runner when
     * {@code agent.runner} is unset, so a container requirement can otherwise read as coming from
     * nowhere. Points at the escape hatch the user is usually reaching for -- a different runner
     * plugin, not a different config key.
     */
    private static String autoSelectedRunnerHint(String runner) {
        return ". Runner `${runner}` was selected automatically as the only agent runner plugin installed - set `agent.runner` to choose another, e.g. `agent.runner = 'langchain4j'` to call the model from the driver process, which needs no container"
    }

    /**
     * Resolve the address a containerized agent task will reach the driver's RPC broker on, and
     * reject the run when the ladder says there is none -- a THIN CALLER of
     * {@link AgentRpcHostResolver}, which owns every row and every rejection, so this guard and the
     * address the endpoint actually advertises are the same decision rather than two that agree by
     * inspection.
     *
     * <p>The rejection stays PRE-IGNITION for the rows that cannot work (E1..E7). That timing is a
     * requirement, not a nicety: a capability is only ever released by its one-hour
     * {@code agent.rpc.capabilityTimeout}, so a configuration allowed to submit would fail only once
     * the job had queued and an instance had booted, holding its request the whole time.
     */
    @PackageScope
    static AgentRpcHost requireBrokerHost(String agentName, String executor, Executor agentExecutor,
            ContainerConfig containerConfig, Object containerOptions, AgentRpcConfig rpc, Session session) {
        final resolved = AgentRpcHostResolver.resolve(agentExecutor, executor, containerConfig, containerOptions, rpc, session)
        if( resolved.resolved )
            return resolved
        throw new ScriptRuntimeException("Agent `${agentName}` on executor `${executor ?: AgentConfig.DEFAULT_EXECUTOR}` cannot determine an address for the driver's agent RPC broker - ${resolved.error}")
    }

    /**
     * The container run options an agent task needs on top of those it declares, for the engine that
     * will run it: Linux Docker (>= 20.10) resolves {@code host.docker.internal} -- the address the
     * task dials the driver's broker on when {@code agent.rpc.remoteHost} is not set -- ONLY when the
     * container is run with {@code --add-host=host.docker.internal:host-gateway}. Docker Desktop,
     * where the name is built in, accepts the same flag, so add it unconditionally for docker rather
     * than leaving the most common local path to fail with a connection timeout. Podman needs
     * nothing: it provides {@code host.containers.internal} itself.
     *
     * <p>What matters is the NAME the task resolves, not whether {@code agent.rpc.remoteHost} was
     * set: spelling the docker alias out explicitly -- as a migrated config, or an overlay
     * inheriting it from a Kubernetes profile, naturally does -- needs the same mapping as leaving
     * it to the default. Only a different address makes the option pointless.
     *
     * <p>A dynamic (closure) or map-valued {@code containerOptions} is resolved per task by
     * {@link nextflow.processor.TaskConfig} and cannot be appended to here, so it is returned
     * unchanged and the flag stays the user's to add.
     */
    static Object withDockerHostGateway(Object containerOptions, String engine, String remoteHost) {
        if( engine != 'docker' )
            return containerOptions
        if( remoteHost && remoteHost != AgentRpcConfig.hostAliasFor('docker') )
            return containerOptions
        if( containerOptions == null )
            return DOCKER_HOST_GATEWAY
        if( containerOptions instanceof CharSequence )
            return "${containerOptions} ${DOCKER_HOST_GATEWAY}".toString()
        log.warn "Agent `containerOptions` is not a plain string, so `${DOCKER_HOST_GATEWAY}` cannot be appended to it - add it yourself, or set `agent.rpc.remoteHost`, if the agent task cannot reach the driver on Linux Docker"
        return containerOptions
    }

    private static final String DOCKER_HOST_GATEWAY = '--add-host=host.docker.internal:host-gateway'

    /**
     * Per-task re-check of the image the canonical launch command needs. {@code agent.container} may
     * be a lazy value (a closure or a GString) that is truthy when {@link #requireContainerized}
     * runs pre-ignition yet resolves to null for the task, so the build-time guard alone cannot
     * keep the in-image paths from being run on the host.
     */
    static void requireTaskContainer(String agentName, String runner, Object container) {
        if( !hasContainer(container) )
            throw new ScriptRuntimeException("Agent `${agentName}` resolved no container image - set `agent.container` to the `${runner}` runner image, which carries the agent proxy and harness the task is launched from")
    }
}
