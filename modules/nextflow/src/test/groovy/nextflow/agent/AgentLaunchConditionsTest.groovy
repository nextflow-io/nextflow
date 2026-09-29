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

import nextflow.Session
import nextflow.SysEnv
import nextflow.agent.rpc.AgentRpcHostResolverTest.TestProbes
import nextflow.container.DockerConfig
import nextflow.container.PodmanConfig
import nextflow.container.SingularityConfig
import nextflow.container.SmolVmConfig
import nextflow.exception.ScriptRuntimeException
import nextflow.executor.AbstractGridExecutor
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import spock.lang.Specification
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.rpc.AgentRpcConfig

class AgentLaunchConditionsTest extends Specification {

    def setup() {
        // AgentRpcConfig reads the environment in its constructor, and a NXF_AGENT_RPC_REMOTE_HOST
        // exported in the developer's shell would silently answer the second rung of the ladder for
        // every feature below, hiding the rung actually under test
        SysEnv.push([:])
    }

    def cleanup() {
        SysEnv.pop()
        AgentRpcHostResolver.reset()
    }

    /**
     * Place the driver: the host facts the address ladder observes, injected for the session the
     * guard will resolve against. These features assert which RUNG the guard takes, so they must not
     * depend on the machine the suite runs on -- a suite running in a container would otherwise see
     * {@code /.dockerenv} and take the containerized-driver row for every alias case below.
     */
    private static Session withDriverHost(Session session, Map probes = [:]) {
        AgentRpcHostResolver.install(session, new TestProbes([outbound: '10.0.3.17', interfaces: ['10.0.3.17']] + probes))
        return session
    }

    /** The `agent.rpc` scope, whose only key any of these features cares about is `remoteHost`. */
    private static AgentRpcConfig rpc(String remoteHost = null) {
        return new AgentRpcConfig(remoteHost != null ? [remoteHost: remoteHost] : [:])
    }

    def 'containerization requires an image AND a container-native executor or an enabled engine'() {
        expect:
        AgentLaunchConditions.willContainerize(container, containerNative, engineEnabled) == expected

        where:
        container | containerNative | engineEnabled || expected
        null      | false           | true          || false
        false     | false           | true          || false
        null      | true            | false         || false
        'image'   | false           | false         || false
        'image'   | false           | true          || true
        'image'   | true            | false         || true
    }

    def 'a canonical launch is rejected whenever the configuration would not containerize it'() {
        when:
        AgentLaunchConditions.requireContainerized('qa', 'pi', executor, container, containerNative, engineEnabled)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('agent.container')
        and: 'the message names the runner image, so it is actionable'
        e.message.contains('`pi` runner image')
        and:
        e.message.contains(reason)

        where:
        executor | container | containerNative | engineEnabled || reason
        // no image at all, whether local or offloaded, and whether or not a container would run
        'local'  | null      | false           | false         || 'must declare a container'
        'local'  | null      | false           | true          || 'must declare a container'
        'k8s'    | false     | true            | false         || 'must declare a container'
        // an image, but nothing that would actually run it in a container: a non-local executor is
        // NOT container-native by itself (the former isOffloaded assumption), and a local executor
        // with no engine enabled would run the in-image paths on the host
        'k8s'    | 'image'   | false           | false         || 'would not run it in a container'
        'slurm'  | 'image'   | false           | false         || 'would not run it in a container'
        'local'  | 'image'   | false           | false         || 'would not run it in a container'
    }

    def 'a rejection says where an unchosen runner came from'() {
        when: 'the user never set `agent.runner`, so the sole installed one was picked for them'
        AgentLaunchConditions.requireContainerized('qa', 'pi', 'local', container, false, false, true)

        then: 'the message names the runner AND accounts for it'
        def err = thrown(ScriptRuntimeException)
        err.message.contains(expected)
        err.message.contains('selected automatically as the only agent runner plugin installed')
        err.message.contains("agent.runner = 'langchain4j'")

        when: 'the runner was chosen explicitly there is nothing to explain'
        AgentLaunchConditions.requireContainerized('qa', 'pi', 'local', container, false, false)

        then:
        def chosen = thrown(ScriptRuntimeException)
        chosen.message.contains(expected)
        !chosen.message.contains('selected automatically')

        where:
        container || expected
        null      || 'must declare a container'
        'image'   || 'would not run it in a container'
    }

    def 'a containerized launch is accepted when the executor or the engine provides the container'() {
        when:
        AgentLaunchConditions.requireContainerized('qa', 'pi', executor, 'image', containerNative, engineEnabled)

        then:
        noExceptionThrown()

        where:
        executor | containerNative | engineEnabled
        'local'  | false           | true
        'local'  | true            | false
        'k8s'    | true            | false
    }

    def 'only the local executor drives a container engine on the driver host'() {
        expect:
        AgentLaunchConditions.isDriverHostEngine(executor) == expected

        where:
        executor | expected
        null     | true
        'local'  | true
        'slurm'  | false
        'k8s'    | false
    }

    def 'the broker host is resolved from the executor, the engine and the driver host itself'() {
        given: 'the executor INSTANCE, which is what the grid row keys off - never a name list'
        def session = withDriverHost(new Session())
        def agentExecutor = grid ? Mock(AbstractGridExecutor) : null

        expect: 'the guard returns the address it resolved, so the caller decides run options from it'
        def host = AgentLaunchConditions.requireBrokerHost('qa', executor, agentExecutor, containerConfig, null, rpc(remoteHost), session)
        host.host == expected
        host.source == source

        where:
        executor | grid  | containerConfig            | remoteHost   || expected                   | source
        // an explicit value always suffices, whatever the executor or engine
        'k8s'    | false | new DockerConfig([:])      | 'driver.svc' || 'driver.svc'               | 'agent.rpc.remoteHost'
        // docker and podman synthesize a name for the host their containers run on
        'local'  | false | new DockerConfig([:])      | null         || 'host.docker.internal'     | 'docker host alias'
        'local'  | false | new PodmanConfig([:])      | null         || 'host.containers.internal' | 'podman host alias'
        // singularity creates no network namespace, so the task IS in the driver's
        'local'  | false | new SingularityConfig([:]) | null         || '127.0.0.1'                | 'host network namespace'
        // and where the task runs elsewhere -- grid, k8s, or an executor with no rung of its
        // own -- the driver's own address on its default route, on one shared row
        'slurm'  | true  | new DockerConfig([:])      | null         || '10.0.3.17'                | 'inferred from default route'
        'k8s'    | false | new DockerConfig([:])      | null         || '10.0.3.17'                | 'inferred from default route'
        'nonesuch'| false| new DockerConfig([:])      | null         || '10.0.3.17'                | 'inferred from default route'
    }

    def 'the broker host is rejected, before ignition, for the configurations no address can serve'() {
        given: 'a driver placed where the ladder has an error row for it'
        def session = withDriverHost(new Session(sessionConfig), probes)

        when:
        AgentLaunchConditions.requireBrokerHost('qa', executor, null, containerConfig, containerOptions, rpc(), session)

        // pre-ignition, because a capability is only ever released by its one-hour timeout: a
        // configuration allowed to submit would fail only once the job had queued, an hour later
        then: 'the run is rejected, and by name'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('`qa`')
        e.message.contains(executor)
        and: 'and the message says which row rejected it, not merely that something did'
        e.message.contains(reason)

        where:
        executor | sessionConfig | containerConfig                    | containerOptions | probes                                  || reason
        // E2: a microVM created with no network at all can dial no address whatsoever
        'local'  | [:]           | new SmolVmConfig([network: false]) | null             | [:]                                     || '`smolvm.network = true`'
        // E7: nothing routable to advertise - a guess here costs an hour of silence
        'k8s'    | [:]           | null                               | null             | [outbound: null, interfaces: []]        || 'no address the agent task could reach the driver on'
    }

    def 'the docker host-gateway run option is added exactly when the docker alias is what the task uses'() {
        expect:
        AgentLaunchConditions.withDockerHostGateway(options, engine, remoteHost) == expected

        where:
        options       | engine        | remoteHost   || expected
        // Linux Docker resolves `host.docker.internal` only with this option; Docker Desktop,
        // where the name is built in, accepts it too - so it is added unconditionally for docker
        null          | 'docker'      | null         || '--add-host=host.docker.internal:host-gateway'
        // an option the agent declares is preserved, never replaced
        '--cpus 2'    | 'docker'      | null         || '--cpus 2 --add-host=host.docker.internal:host-gateway'
        // spelling the docker alias out explicitly needs the mapping just as much as the default
        // does - a migrated config, or an overlay inheriting the value from a k8s profile, does
        null          | 'docker'      | 'host.docker.internal' || '--add-host=host.docker.internal:host-gateway'
        // nothing to add: the task was given an explicit address for the driver
        null          | 'docker'      | 'driver.svc' || null
        '--cpus 2'    | 'docker'      | 'driver.svc' || '--cpus 2'
        // podman provides `host.containers.internal` itself, and no other engine takes the option
        null          | 'podman'      | null         || null
        null          | 'singularity' | null         || null
        // a container-native executor: no local engine runs the container at all
        '--cpus 2'    | null          | 'driver.svc' || '--cpus 2'
    }

    def 'a dynamic containerOptions value is left untouched rather than stringified'() {
        given: 'a closure is resolved per task by TaskConfig, so it cannot be appended to here'
        final dynamic = { -> '--cpus 2' }

        expect:
        AgentLaunchConditions.withDockerHostGateway(dynamic, 'docker', null).is(dynamic)
    }

    /**
     * A session whose executor answers the three questions the guard asks it -- does it manage
     * containers itself, which engine block must be enabled, and (through its TYPE) where it runs
     * the task. Stubbing the EXECUTOR rather than passing the answers in is what puts the resolution
     * under test: which executor is interrogated, and which engine block is read for it.
     *
     * <p>The driver host is placed too, so the address rung a feature exercises is the one it asked
     * for and not one the machine running the suite happens to supply.
     */
    private Session sessionWith(Map config, boolean containerNative, String configEngine,
            Class<? extends Executor> executorType = Executor, Map probes = [:]) {
        final executor = Stub(executorType) {
            isContainerNative() >> containerNative
            containerConfigEngine() >> configEngine
        }
        final session = new Session(config)
        session.executorFactory = Stub(ExecutorFactory) {
            getExecutorByName(_, _) >> executor
        }
        return withDriverHost(session, probes)
    }

    def 'an agent with no container is rejected before the executor is even resolved'() {
        given: 'a factory that blows up if it is asked for an executor at all'
        def session = new Session([docker: [enabled: true]])
        session.executorFactory = Stub(ExecutorFactory) {
            getExecutorByName(_, _) >> { throw new IllegalStateException('the executor must not be resolved') }
        }

        when: 'no image is declared - not even an executor that cannot be instantiated hides that'
        AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'local', container, null, rpc('driver.internal'), session)

        then: 'the message names the directive to set AND what the runner image is for'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare a container')
        e.message.contains('agent.container')
        e.message.contains('`pi` runner image')

        where: 'both an absent image and the explicit opt-out'
        container << [ null, false ]
    }

    def 'an agent with a container is rejected when nothing would run it in one'() {
        given: 'a grid executor: not container-native, and no engine enabled in the config'
        def session = sessionWith([:], false, null)

        when:
        AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'slurm', 'agent-image:test', null, rpc('driver.internal'), session)

        then: 'the in-image paths would be exec-ed on the host, so fail here instead'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('would not run it in a container')
        e.message.contains('agent.container')
        e.message.contains('`pi` runner image')
    }

    def 'a container-native executor containerizes the agent with no engine enabled'() {
        given: 'k8s: container-native, and it reads the docker config block'
        def session = sessionWith([:], true, 'docker')

        when: 'no `docker.enabled` anywhere - the executor manages containers itself'
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'k8s', 'agent-image:test', null, rpc('driver.svc'), session)

        then:
        noExceptionThrown()
        and: 'no engine is reported: the container is launched in the cluster, not by the driver'
        launch.containerEngine == null
        launch.brokerHost.host == 'driver.svc'

        when: 'and with no address given, the driver advertises the address it routes out on'
        def inferred = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'k8s', 'agent-image:test', null, rpc(), session)

        then: 'the docker host alias names nothing a pod can resolve, so it must not stand in'
        inferred.brokerHost.host == '10.0.3.17'
        inferred.brokerHost.source == 'inferred from default route'
    }

    def 'a local executor with an enabled engine containerizes the agent'() {
        given: 'the local executor is not container-native and pins no engine of its own'
        def session = sessionWith([docker: [enabled: true]], false, null)

        when: 'no explicit `agent.rpc.remoteHost` - docker names the container host itself'
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'local', 'agent-image:test', null, rpc(), session)

        then:
        noExceptionThrown()
        and: 'the engine is reported back, because it runs on the driver host and takes run options'
        launch.containerEngine == 'docker'
        and: 'together with the address it resolved, which is what decides those run options'
        launch.brokerHost.host == 'host.docker.internal'
        launch.brokerHost.source == 'docker host alias'
    }

    def 'a Fusion-enabled local executor still reaches the driver at the docker host alias'() {
        given: 'Fusion makes the LOCAL executor container-native, but docker still runs on the driver host'
        def session = sessionWith([docker: [enabled: true]], true, null)

        when: 'no explicit `agent.rpc.remoteHost`'
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'local', 'agent-image:test', null, rpc(), session)

        then: 'container-native does not mean remote - the alias applies and no address is demanded'
        noExceptionThrown()
        launch.containerEngine == 'docker'
        launch.brokerHost.host == 'host.docker.internal'
    }

    def 'a remote executor driving an engine of its own reaches the driver at its own address'() {
        given: 'a grid executor with docker enabled: the daemon runs on the compute node, not the driver'
        def session = sessionWith([docker: [enabled: true]], false, null, AbstractGridExecutor)

        when:
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'slurm', 'agent-image:test', null, rpc(), session)

        // the node's `host.docker.internal` names that node's host, never the driver, so what the
        // job dials back on is the driver's own address on its default route
        then: 'the address is inferred rather than demanded'
        launch.brokerHost.host == '10.0.3.17'
        launch.brokerHost.source == 'inferred from default route'
        and: 'and no engine is reported: the driver adds no run options to a remote launch'
        launch.containerEngine == null

        when: 'the driver address is declared, it still wins over what was inferred'
        def explicit = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'slurm', 'agent-image:test', null, rpc('driver.internal'), session)

        then:
        explicit.brokerHost.host == 'driver.internal'
        explicit.containerEngine == null
    }

    def 'the enabled engine is discovered through the executor and names the container host'() {
        given: 'podman rather than docker, discovered because the local executor pins no engine'
        def session = sessionWith([podman: [enabled: true]], false, null)

        when: 'no explicit `agent.rpc.remoteHost` - podman has a host alias of its own'
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'local', 'agent-image:test', null, rpc(), session)

        then:
        noExceptionThrown()
        launch.containerEngine == 'podman'
        launch.brokerHost.host == 'host.containers.internal'
    }

    def 'the engine block consulted is the one the EXECUTOR asks for, not whichever is enabled'() {
        given: 'an executor that pins the docker block (as k8s does) while podman is the enabled engine'
        def session = sessionWith([podman: [enabled: true]], false, 'docker')

        when:
        AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'custom', 'agent-image:test', null, rpc('driver.internal'), session)

        then: 'an engine the executor does not use cannot satisfy the guard - the task would not see it'
        // the engine-agnostic session.getContainerConfig() would accept this, which is exactly the
        // disagreement with the task's own view (TaskRun.getContainerConfig asks the executor too)
        def e = thrown(ScriptRuntimeException)
        e.message.contains('would not run it in a container')
    }

    def 'a containerizing engine with no host alias reaches the driver in its own network namespace'() {
        given: 'singularity DOES containerize the task, it just has no name for the container host'
        def session = sessionWith([singularity: [enabled: true]], false, null)

        when: 'no address is given: singularity creates no network namespace of its own'
        def launch = AgentLaunchConditions.requireCanonicalLaunch('qa', 'pi', 'local', 'agent-image:test', null, rpc(), session)

        then: 'the task is IN the driver`s namespace, so the driver`s loopback is what it dials'
        launch.brokerHost.host == '127.0.0.1'
        launch.brokerHost.source == 'host network namespace'
        launch.containerEngine == 'singularity'

    }
}
