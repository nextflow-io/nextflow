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

import nextflow.agent.AgentLaunchSpec
import nextflow.agent.AgentRunnerRequest
import spock.lang.Specification

/**
 * Pins the LIVE production path: {@code getLaunchSpec()} is what an agent task actually execs.
 */
class PiAgentLaunchSpecTest extends Specification {

    def 'the launch spec describes the in-image command only'() {
        given:
        def spec = new PiAgentRunner().getLaunchSpec()

        expect: 'the proxy, then the packaged harness after the separator, at their IN-IMAGE paths'
        // these must match the paths the runner image is built with (see the plugin Dockerfile);
        // a canonical agent task always runs in that image, so there is no driver-local variant
        spec.command() == ['/usr/local/bin/agent-rpc', '--', 'node', '/opt/nf-agent-pi/runner.mjs']
    }

    def 'the runner has no driver-local execution path'() {
        when: 'the deleted stdio route is invoked'
        new PiAgentRunner().run(new AgentRunnerRequest(prompt: 'Q'))

        then: 'it fails loudly instead of half-working, since the runtime lives in the image'
        def e = thrown(UnsupportedOperationException)
        e.message.contains('container task')
    }

    def 'proxy arguments are composed without parsing the command separator'() {
        given:
        def spec = new AgentLaunchSpec(
            ['launcher', '--', 'agent-rpc'],
            ['node', '/runner.mjs'])

        expect: 'an existing proxy argument named `--` is not mistaken for the harness separator'
        spec.command(['--endpoint', 'driver:1234']) == [
            'launcher', '--', 'agent-rpc',
            '--endpoint', 'driver:1234',
            '--',
            'node', '/runner.mjs' ]
    }
}
