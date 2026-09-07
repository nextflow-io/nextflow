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

import spock.lang.Specification

class AgentProtocolSpecTest extends Specification {

    def 'should create the portable protocol payload from a runner request'() {
        given:
        final request = new AgentRunnerRequest(
            model: 'provider/model',
            instruction: 'instruction',
            goal: 'goal',
            prompt: 'prompt',
            inputJson: '{"value":1}',
            outputSchema: [type: 'string'],
            toolSpecs: [],
            nativeToolNames: ['read'],
            skills: [],
            maxIterations: 0,
            trace: true,
            temperature: null,
            workDir: '.')

        expect:
        AgentProtocolSpec.fromRequest(request) == [
            model: 'provider/model',
            instruction: 'instruction',
            goal: 'goal',
            prompt: 'prompt',
            inputJson: '{"value":1}',
            outputSchema: [type: 'string'],
            toolSpecs: [],
            nativeToolNames: ['read'],
            skills: [],
            maxIterations: 20,
            trace: true,
            temperature: null,
            workDir: '.',
            baseUrl: null ]
    }

    // -----------------------------------------------------------------------
    // The runner split. `toolSpecs` carries the BROKERED tools; the runner-native ones travel
    // beside it as bare names the runner enables from its own tool set.
    // -----------------------------------------------------------------------

    def 'should carry the runner-native names beside the brokered descriptors, never inside them'() {
        given:
        final request = new AgentRunnerRequest(
            model: 'openai/gpt-4o',
            prompt: 'p',
            toolSpecs: [new ToolDescriptor('SAMTOOLS_SORT', 'sort it', [type: 'object'], null)],
            nativeToolNames: ['read', 'write', 'bash'])

        when:
        final payload = AgentProtocolSpec.fromRequest(request)

        then: 'a native tool has no descriptor: it is a name the runner recognises as its own'
        payload.toolSpecs*.name == ['SAMTOOLS_SORT']
        payload.nativeToolNames == ['read', 'write', 'bash']
    }

    def 'should refuse to build a payload whose two tool halves overlap'() {
        given: 'a name claimed by BOTH halves -- the runner would be told to serve it itself AND'
        // to call the driver back for it, and the broker would authorize the callback
        final request = new AgentRunnerRequest(
            model: 'openai/gpt-4o',
            prompt: 'p',
            toolSpecs: [new ToolDescriptor('read', 'a process named read', [type: 'object'], null)],
            nativeToolNames: ['read', 'bash'])

        when:
        AgentProtocolSpec.fromRequest(request)

        then: 'it fails HERE, before any frame is written, rather than in the container'
        final err = thrown(IllegalStateException)
        err.message.contains('partition violated')
        err.message.contains('read')
        !err.message.contains('bash')
    }

    // -----------------------------------------------------------------------
    // Credential containment (design D6). This payload crosses the plaintext gRPC link to a
    // possibly remote agent task, so the endpoint travels and the credential must not.
    // -----------------------------------------------------------------------

    def 'should carry the resolved endpoint but NEVER the credential'() {
        given: 'a request as core builds it: both the endpoint and the credential resolved'
        final request = new AgentRunnerRequest(
            model: 'openai/llama-3.3-70b',
            prompt: 'prompt',
            apiKey: 'sk-must-not-travel-9f2c',
            baseUrl: 'http://localhost:8000/v1')

        when:
        final payload = AgentProtocolSpec.fromRequest(request)

        then: 'the endpoint travels -- a remote runner must target the endpoint the DRIVER resolved'
        payload.baseUrl == 'http://localhost:8000/v1'

        and: 'the credential does not, under any key'
        !payload.containsKey('apiKey')
        !payload.values().contains('sk-must-not-travel-9f2c')
        !payload.toString().contains('sk-must-not-travel-9f2c')
    }

    def 'the payload key set carries nothing credential-shaped, whatever the transport adds beside it'() {
        given: 'design D4 now DOES send the credential to a pi task -- as a top-level start-frame'
        // field the broker adds under TLS, never inside this map. The invariant therefore has to be
        // asserted on the SHAPE, not merely on the one key: this payload is the half a transport is
        // free to relay verbatim, log or persist, and it must stay credential-free by construction.
        final request = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            instruction: 'i', goal: 'g', prompt: 'p',
            inputJson: '{}', outputSchema: [type: 'string'],
            toolSpecs: [], skills: [], maxIterations: 3, trace: true,
            temperature: 0.0d, workDir: '/w',
            apiKey: 'sk-must-not-travel-1c4b',
            baseUrl: 'https://api.openai.com/v1')

        when:
        final payload = AgentProtocolSpec.fromRequest(request)

        then: 'no key name suggests a credential, so a later field cannot slip one in unnoticed'
        payload.keySet().every { !(it ==~ /(?i).*(apikey|api_key|secret|token|password|credential).*/) }
        and: 'and no VALUE is the credential, under any key'
        !payload.values().contains('sk-must-not-travel-1c4b')

        and: 'the class declares no credential field either'
        AgentProtocolSpec.declaredFields.every { !(it.name ==~ /(?i).*(apikey|api_key|secret|token|credential).*/) }
    }

    def 'the credential leaks into neither the request toString nor the persisted task info'() {
        given:
        final request = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            prompt: 'p',
            apiKey: 'sk-must-not-log-4a71',
            baseUrl: 'http://gateway/v1')

        expect: '@ToString(excludes=apiKey): an interpolated request cannot leak the key into .nextflow.log'
        !request.toString().contains('sk-must-not-log-4a71')
        and: 'the endpoint is not a secret and stays visible for diagnostics'
        request.toString().contains('http://gateway/v1')

        and: 'AgentTaskInfo -- persisted by the lineage observer -- declares no credential field at all'
        AgentTaskInfo.declaredFields.every { !(it.name ==~ /(?i).*(apikey|api_key|secret|token|password).*/) }
    }
}
