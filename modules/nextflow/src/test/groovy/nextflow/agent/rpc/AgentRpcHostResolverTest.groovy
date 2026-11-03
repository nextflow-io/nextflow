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

import nextflow.SysEnv
import nextflow.container.AppleContainerConfig
import nextflow.container.ApptainerConfig
import nextflow.container.CharliecloudConfig
import nextflow.container.ContainerConfig
import nextflow.container.DockerConfig
import nextflow.container.PodmanConfig
import nextflow.container.SarusConfig
import nextflow.container.ShifterConfig
import nextflow.container.SingularityConfig
import nextflow.container.SmolVmConfig
import nextflow.executor.AbstractGridExecutor
import nextflow.executor.Executor
import spock.lang.Specification

/**
 * Every row of the host-resolution ladder, and every error row, driven through the injected probes
 * so no test opens a socket, reads {@code /proc}, or needs a container, a cluster or a cloud
 * instance. The probes are the ONLY way the resolver observes its host, which is what makes this
 * possible -- a row that could only be reached by really being on that host would not be testable.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentRpcHostResolverTest extends Specification {

    /**
     * A recording stand-in for the host facts. Every answer is settable and every call is counted,
     * so a test can both place the driver in a situation and assert the resolver did not ask twice.
     */
    static class TestProbes implements Probes {
        String outbound = '10.0.3.17'
        List<String> interfaces = ['10.0.3.17']

        final Map<String,Integer> calls = [:].withDefault { 0 }

        private def count(String name) { calls[name] = calls[name] + 1 }

        @Override String outboundAddress() { count('outboundAddress'); outbound }
        @Override List<String> interfaceAddresses() { count('interfaceAddresses'); interfaces }
    }

    def setup() {
        // AgentRpcConfig reads the environment in its constructor, and NXF_AGENT_RPC_REMOTE_HOST
        // exported in the developer's shell would silently answer R2 for every row below
        SysEnv.push([:])
    }

    def cleanup() {
        SysEnv.pop()
        AgentRpcHostResolver.reset()
    }

    private AgentRpcHostResolver resolver(TestProbes probes, Map sessionConfig = [:]) {
        return new AgentRpcHostResolver(sessionConfig, probes)
    }

    private static AgentRpcConfig rpc(Map opts = [:]) {
        return new AgentRpcConfig(opts)
    }

    /** The local executor with a docker engine, the shape every R3/R4/R5 row starts from. */
    private static ContainerConfig docker(Map opts = [:]) { new DockerConfig(opts) }

    // =======================================================================================
    // R1, R2 -- the explicit rungs
    // =======================================================================================

    def 'R1: an explicit `agent.rpc.remoteHost` wins over every inferred row, and probes nothing'() {
        given:
        def probes = new TestProbes(outbound: '172.17.0.5')

        when: 'a configuration that would otherwise resolve R3'
        def host = resolver(probes).resolve(null, 'local', docker(), null, rpc([remoteHost: 'driver.internal']))

        then:
        host.resolved
        host.host == 'driver.internal'
        host.source == 'agent.rpc.remoteHost'

        and: 'no probe ran at all: an explicit value is not something to second-guess'
        probes.calls.isEmpty()
    }

    def 'R2: the environment answers when the config does not, and says so'() {
        given:
        def probes = new TestProbes()
        SysEnv.push([NXF_AGENT_RPC_REMOTE_HOST: 'driver.from.env'])

        when:
        def host = resolver(probes).resolve(null, 'k8s', null, null, rpc())

        then:
        host.host == 'driver.from.env'
        host.source == 'NXF_AGENT_RPC_REMOTE_HOST'
        probes.calls.isEmpty()

        cleanup:
        SysEnv.pop()
    }

    // =======================================================================================
    // R3 -- the containerized driver (a fix to shipped behaviour, not a new capability)
    // =======================================================================================

    // =======================================================================================
    // R4 -- the container shares the driver's network namespace
    // =======================================================================================

    def 'R4: the engines that create no network namespace resolve to loopback'() {
        expect:
        def host = resolver(new TestProbes()).resolve(null, 'local', config, null, rpc())
        host.host == '127.0.0.1'
        host.source == 'host network namespace'

        where:
        config << [
            new SingularityConfig([:]),
            new ApptainerConfig([:]),
            // an explicit `--network host` is not a namespace of its own
            new ApptainerConfig([runOptions: '--network host']),
        ]
    }

    // =======================================================================================
    // R5, R6 -- the local engine rows
    // =======================================================================================

    def 'R5: a non-containerized driver on a local daemon gets the engine host alias'() {
        expect:
        def host = resolver(new TestProbes()).resolve(null, 'local', config, null, rpc())
        host.host == expected
        host.source == source

        where:
        config                | expected                      | source
        docker()              | 'host.docker.internal'        | 'docker host alias'
        new PodmanConfig([:]) | 'host.containers.internal'    | 'podman host alias'
    }

    def 'R6: an isolated engine with no alias is given the default-route address'() {
        given:
        def probes = new TestProbes(outbound: '192.168.1.131', interfaces: ['192.168.1.131'])

        expect:
        def host = resolver(probes).resolve(null, 'local', config, null, rpc())
        host.host == '192.168.1.131'
        host.source == 'inferred from default route'

        where:
        config << [new AppleContainerConfig([:]), new SmolVmConfig([:]), new SmolVmConfig([network: true])]
    }

    // =======================================================================================
    // R7, R8, R9 -- the executor rows
    // =======================================================================================

    def 'R7: a grid executor is recognised by its TYPE, not by a name list'() {
        given: 'a subclass whose @ServiceName no hand-written list would carry'
        def executor = Mock(AbstractGridExecutor)
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3'])

        when:
        def host = resolver(probes).resolve(executor, 'nqsii', docker(), null, rpc())

        then:
        host.host == '10.1.2.3'
        host.source == 'inferred from default route'
    }

    def 'R7: a multi-homed submit node still resolves, but says another address exists'() {
        given: 'the campus NIC owns the default route; the compute fabric is the other one'
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3', '192.168.40.7'])

        when:
        def host = resolver(probes).resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then:
        host.host == '10.1.2.3'
        host.warnings.size() == 1
        host.warnings[0].contains('192.168.40.7')
        host.warnings[0].contains('agent.rpc.remoteHost')
    }

    def 'R7: a k8s driver is given its own address, on the same rung as every other remote executor'() {
        given: 'a driver in the cluster the agent pods run in'
        def probes = new TestProbes(outbound: '10.42.0.9', interfaces: ['10.42.0.9'])

        when:
        def host = resolver(probes).resolve(Mock(Executor), 'k8s', null, null, rpc())

        then:
        host.host == '10.42.0.9'
        host.source == 'inferred from default route'
        and: 'a hostNetwork pod yields the NODE address, which is correct but is not a pod IP'
        !host.source.contains('pod IP')
    }

    // =======================================================================================
    // error rows -- one message each, naming what was tried and what to set
    // =======================================================================================

    def 'E2: a smolvm microVM with no network at all is rejected'() {
        when:
        def host = resolver(new TestProbes()).resolve(null, 'local', new SmolVmConfig([network: false]), null, rpc())

        then:
        !host.resolved
        host.code == 'E2'
        host.error.contains('smolvm.network')
    }

    def 'E7: an address that is not routable from another host is rejected, naming every step tried'() {
        given:
        def probes = new TestProbes(outbound: address, interfaces: [])

        when:
        def host = resolver(probes).resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then:
        !host.resolved
        host.code == 'E7'
        and: 'the message is a trace of the ladder, not a bare "cannot determine"'
        host.error.contains('agent.rpc.remoteHost')
        host.error.contains('NXF_AGENT_RPC_REMOTE_HOST')
        host.error.contains('default route')

        where:
        address << ['127.0.0.1', '0.0.0.0', '169.254.12.9', null]
    }

    def 'a host with no default route falls back to its ONE routable interface address'() {
        given: 'an air-gapped submit node: one NIC on the compute fabric, no default route'
        def probes = new TestProbes(outbound: null, interfaces: ['10.20.0.5'])

        when:
        def host = resolver(probes).resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then: 'the address the compute nodes reach it on was enumerated all along'
        host.resolved
        host.host == '10.20.0.5'
        host.source.contains('no default route')
    }

    def 'a host with no default route falls back to its ONE PUBLIC address, and says so'() {
        given: 'no default route, and a private plus a public address on the interfaces'
        def probes = new TestProbes(outbound: null, interfaces: ['10.0.3.17', '203.0.113.9'])

        when: 'a remote executor, i.e. the deployment where a public address is the fair fallback'
        def host = resolver(probes).resolve(Mock(Executor), 'k8s', null, null, rpc())

        then: 'the public one settles an otherwise ambiguous set'
        host.resolved
        host.host == '203.0.113.9'

        and: 'and the operator is told, because a public address is a consequential choice'
        host.warnings.any { it.contains('public address') }
    }

    def 'a public fallback is NOT taken when two public addresses are equally plausible'() {
        given:
        def probes = new TestProbes(outbound: null, interfaces: ['203.0.113.9', '198.51.100.4'])

        when:
        def host = resolver(probes).resolve(Mock(Executor), 'k8s', null, null, rpc())

        then: 'advertising a guess is the failure this design exists to avoid'
        !host.resolved
        host.code == 'E7'
    }

    def 'a host with no default route and SEVERAL routable addresses is rejected, not guessed at'() {
        given:
        def probes = new TestProbes(outbound: null, interfaces: ['10.20.0.5', '192.168.9.4'])

        when:
        def host = resolver(probes).resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then: 'nothing says which one the agent task can route to, and a guess is the failure mode'
        !host.resolved
        host.code == 'E7'
        host.error.contains('10.20.0.5')
        host.error.contains('192.168.9.4')
        host.error.contains('agent.rpc.remoteHost')
    }

    def 'the multi-homed warning ignores the SAME NIC in another address family'() {
        given: 'the ordinary dual-stack Linux host: one NIC, an IPv4 and a global IPv6 address'
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3', '2a01:4f8:1:2::1'])

        when:
        def host = resolver(probes).resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then: 'warning on every ordinary driver would train operators to ignore the one that matters'
        host.host == '10.1.2.3'
        host.warnings.isEmpty()
    }

    def 'R7: an executor with no rung of its own takes the same remote row as k8s and grid'() {
        given: 'an engine that names a driver-host address, but an executor that is not on the driver host'
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3'])

        when:
        def host = resolver(probes).resolve(Mock(Executor), 'nonesuch', docker(), null, rpc())

        then: 'the engine alias is NOT used - it names the wrong machine - and the routable address is'
        host.resolved
        host.host == '10.1.2.3'
        host.source == 'inferred from default route'
    }

    def 'R7: a non-namespace engine on the driver host with no rung of its own also falls through'() {
        given: 'shifter/charliecloud/sarus get no special handling; the outbound address serves them'
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3'])

        expect:
        def host = resolver(probes).resolve(null, 'local', config, null, rpc())
        host.resolved
        host.host == '10.1.2.3'

        where:
        config << [new ShifterConfig([:]), new CharliecloudConfig([:]), new SarusConfig([:])]
    }

    // =======================================================================================
    // memoization, and the broker hand-off
    // =======================================================================================

    def 'every probe runs at most once, however many agent definitions resolve'() {
        given:
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3'])
        def resolver = resolver(probes)

        when: 'three agent definitions resolve through the same session resolver'
        3.times { resolver.resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc()) }

        then: 'these answer questions about the HOST, so asking again could only cost'
        probes.calls['outboundAddress'] == 1
        probes.calls['interfaceAddresses'] == 1
    }

    def 'each agent definition keeps its OWN address, so neither can be advertised to the other'() {
        given: 'one agent runs locally on docker, another on a grid executor'
        def probes = new TestProbes(outbound: '10.1.2.3', interfaces: ['10.1.2.3'])
        def resolver = resolver(probes)

        when:
        def local = resolver.resolve(null, 'local', docker(), null, rpc())
        def grid = resolver.resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc())

        then: 'the resolver is stateless across definitions -- the host rides on the request'
        local.host == 'host.docker.internal'
        grid.host == '10.1.2.3'

        and: 'and the order the definitions are built in cannot change either answer'
        resolver.resolve(Mock(AbstractGridExecutor), 'slurm', docker(), null, rpc()).host == '10.1.2.3'
        resolver.resolve(null, 'local', docker(), null, rpc()).host == 'host.docker.internal'
    }

    // =======================================================================================
    // parity with the thin callers
    // =======================================================================================

    def 'an address is usable only when another host could route to it'() {
        expect:
        AgentRpcHostResolver.usableAddress(address) == expected

        where:
        address           | expected
        '10.0.3.17'       | true
        '192.168.1.131'   | true
        '127.0.0.1'       | false      // the driver's own loopback is not the task's
        '0.0.0.0'         | false      // the bind address, never an address to dial
        '169.254.169.254' | false      // link-local: reachable only on the same link
        ''                | false
        null              | false
    }

    def 'AgentRpcConfig.resolveRemoteHost goes through the same ladder'() {
        expect: 'the local rows, which are all an engine name alone can answer'
        rpc().resolveRemoteHost('docker') == 'host.docker.internal'
        rpc().resolveRemoteHost('podman') == 'host.containers.internal'
        rpc().resolveRemoteHost('singularity') == '127.0.0.1'
        rpc().resolveRemoteHost('apptainer') == '127.0.0.1'

        and: 'an explicit value still outranks every row'
        rpc([remoteHost: 'driver.internal']).resolveRemoteHost('docker') == 'driver.internal'
    }
}
