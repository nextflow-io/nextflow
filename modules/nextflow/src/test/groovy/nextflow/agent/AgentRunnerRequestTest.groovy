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

import nextflow.exception.AbortOperationException
import spock.lang.Specification

class AgentRunnerRequestTest extends Specification {

    def 'should build via named args with goal as the last field'() {
        when:
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            instruction: 'sys',
            prompt: 'p',
            maxIterations: 7,
            tools: [],
            outputSchema: null,
            inputJson: '{}',
            toolSpecs: null,
            dispatch: null,
            requestTimeoutSeconds: 30,
            goal: 'reach the objective',
            skills: [new SkillDescriptor('greet', 'a greeting skill', 'say hi', [])])

        then:
        req.model == 'openai/gpt-5-mini'
        req.instruction == 'sys'
        req.prompt == 'p'
        req.maxIterations == 7
        req.inputJson == '{}'
        req.requestTimeoutSeconds == 30
        req.goal == 'reach the objective'
        req.tools == []
        req.outputSchema == null
        req.toolSpecs == null
        req.skills.size() == 1
        req.skills[0].name == 'greet'
    }

    def 'should default skills to null when omitted'() {
        when:
        def req = new AgentRunnerRequest(model: 'm', prompt: 'p')
        then:
        req.skills == null
    }

    def 'should default goal to null when omitted'() {
        when:
        def req = new AgentRunnerRequest(model: 'm', prompt: 'p')
        then:
        req.goal == null
    }

    // -----------------------------------------------------------------------
    // The brokered/runner-native partition. Two fields, disjoint by construction: `toolSpecs` is
    // what the runner may call the DRIVER back for, `nativeToolNames` what it serves itself.
    // -----------------------------------------------------------------------

    def 'should authorize the brokered names only'() {
        given:
        def req = new AgentRunnerRequest(
            model: 'm', prompt: 'p',
            toolSpecs: [
                new ToolDescriptor('GREET', 'greet', [type: 'object'], null),
                new ToolDescriptor('SHOUT', 'shout', [type: 'object'], null) ],
            nativeToolNames: ['read', 'grep', 'bash'])

        expect: 'a runner-native tool never becomes callable over the broker'
        req.brokeredToolNames() == ['GREET', 'SHOUT'] as Set
    }

    def 'should authorize nothing when no tool is declared'() {
        expect:
        new AgentRunnerRequest(model: 'm', prompt: 'p').brokeredToolNames() == [] as Set
    }

    def 'should reject a request claiming one name for both halves'() {
        given: 'a process named `find`, alongside the `fs:find` the runner serves'
        def req = new AgentRunnerRequest(
            model: 'm', prompt: 'p',
            toolSpecs: [new ToolDescriptor('find', 'a process', [type: 'object'], null)],
            nativeToolNames: ['read', 'find'])

        when:
        req.brokeredToolNames()

        then: 'the allowlist is never built: it would relocate a container-side tool into the driver'
        def err = thrown(IllegalStateException)
        err.message.contains('partition violated')
        err.message.contains('find')
    }

    // -----------------------------------------------------------------------
    // The credential a runner presents (design D8). ONE rule, here, because both runners consume
    // it: langchain4j through `credentialFor`, pi through `credential()` on the start frame.
    // -----------------------------------------------------------------------

    def 'credential returns the resolved key when one resolved'() {
        expect:
        new AgentRunnerRequest(model: 'openai/gpt-5-mini', apiKey: 'sk-real').credential() == 'sk-real'
        and: 'an endpoint does not displace it'
        new AgentRunnerRequest(model: 'openai/gpt-5-mini', apiKey: 'sk-real', baseUrl: 'http://local/v1').credential() == 'sk-real'
        and: 'nor does a non-openai provider: the key was scoped to it by core before it got here'
        new AgentRunnerRequest(model: 'anthropic/claude-sonnet-4', apiKey: 'sk-real').credential() == 'sk-real'
    }

    def 'credential falls back to the placeholder for an openai endpoint with no key'() {
        given: 'a local vLLM/Ollama needs no credential, but both runners require SOMETHING:'
        // langchain4j's OpenAI client rejects an empty key, and pi fails the run with
        // `No API key found for openai` -- exactly the local-first case D8 exists to unblock
        expect:
        new AgentRunnerRequest(model: 'openai/llama-3.3-70b', baseUrl: 'http://localhost:8000/v1').credential() ==
            AgentRunnerRequest.PLACEHOLDER_API_KEY

        and: 'with no endpoint there is nothing to talk to, so no placeholder is invented'
        new AgentRunnerRequest(model: 'openai/gpt-5-mini').credential() == null
    }

    def 'the placeholder is never sent for a non-openai provider'() {
        given: 'a runner installs what it is given as the credential OF THAT PROVIDER, and that'
        // ownership beats the ambient environment (pi: setRuntimeApiKey shadows ANTHROPIC_API_KEY),
        // so a placeholder here would MASK a credential the runner can resolve by itself
        expect:
        new AgentRunnerRequest(model: 'anthropic/claude-sonnet-4', baseUrl: 'http://gateway/v1').credential() == null
        and: 'the same holds for a model id with no provider prefix'
        new AgentRunnerRequest(model: 'gpt-5-mini', baseUrl: 'http://gateway/v1').credential() == null
    }

    def 'credentialFor is the same rule for a caller holding only the pair'() {
        expect: 'used by the langchain4j ChatModelFactory, which has already checked the protocol'
        AgentRunnerRequest.credentialFor('sk-real', null) == 'sk-real'
        AgentRunnerRequest.credentialFor('sk-real', 'http://local/v1') == 'sk-real'
        AgentRunnerRequest.credentialFor(null, 'http://local/v1') == AgentRunnerRequest.PLACEHOLDER_API_KEY
        AgentRunnerRequest.credentialFor(null, null) == null

        and: 'the placeholder is not a credential-shaped string: it can never be mistaken for one'
        !AgentRunnerRequest.PLACEHOLDER_API_KEY.startsWith('sk-')
    }

    // -----------------------------------------------------------------------
    // Design D5: the placeholder assumes the endpoint needs no credential, which is true of a
    // local vLLM/Ollama and plainly false of a provider's own API. `api.anthropic.com` is not a
    // local server, so that combination is diagnosed here instead of buying an opaque 401.
    // -----------------------------------------------------------------------

    def 'a well-known provider endpoint with no credential is an error, not a placeholder'() {
        when:
        AgentRunnerRequest.credentialFor(null, ENDPOINT)

        then:
        def error = thrown(AbortOperationException)
        and: 'the message names the provider, the endpoint and the variables the ladder consults'
        error.message.contains(PROVIDER)
        error.message.contains(ENDPOINT)
        error.message.contains('`agent.apiKey`')
        error.message.contains('NXF_AGENT_API_KEY')
        error.message.contains(VAR)

        where:
        ENDPOINT                        | PROVIDER     | VAR
        'https://api.openai.com/v1'     | 'openai'     | 'OPENAI_API_KEY'
        'https://api.anthropic.com/v1'  | 'anthropic'  | 'ANTHROPIC_API_KEY'
        'https://openrouter.ai/api/v1'  | 'openrouter' | 'OPENROUTER_API_KEY'
        'https://api.mistral.ai/v1'     | 'mistral'    | 'MISTRAL_API_KEY'
    }

    def 'the D5 check keys off the HOST, so an unrecognized endpoint keeps the placeholder'() {
        expect: 'a local server, a corporate gateway, a self-hosted mirror -- all still unblocked'
        AgentRunnerRequest.credentialFor(null, 'http://localhost:8000/v1') == AgentRunnerRequest.PLACEHOLDER_API_KEY
        AgentRunnerRequest.credentialFor(null, 'https://gateway.corp/v1') == AgentRunnerRequest.PLACEHOLDER_API_KEY

        and: 'and a provider name that is only in the PATH is not a provider endpoint'
        AgentRunnerRequest.credentialFor(null, 'https://evil.example.com/openai/v1') == AgentRunnerRequest.PLACEHOLDER_API_KEY
        AgentRunnerRequest.credentialFor(null, 'https://api.openai.com.evil.example/v1') == AgentRunnerRequest.PLACEHOLDER_API_KEY

        and: 'a resolved credential short-circuits the check entirely -- there is nothing missing'
        AgentRunnerRequest.credentialFor('sk-real', 'https://api.anthropic.com/v1') == 'sk-real'
    }

    def 'credential() applies D5 only where the placeholder could have been substituted'() {
        when: 'an openai-protocol model points at a provider API with nothing resolved'
        new AgentRunnerRequest(model: 'openai/gpt-5-mini', baseUrl: 'https://api.openai.com/v1').credential()

        then:
        thrown(AbortOperationException)

        when: 'the same endpoint, but the model is not openai-protocol'
        // credential() never reaches credentialFor there: the runner may resolve the key itself
        // (pi reads its own store and the provider variables), so this is not core's call to abort
        def other = new AgentRunnerRequest(model: 'anthropic/claude-sonnet-4', baseUrl: 'https://api.anthropic.com/v1').credential()

        then:
        other == null

        when: 'a credential DID resolve'
        def resolved = new AgentRunnerRequest(model: 'openai/gpt-4o', apiKey: 'sk-real', baseUrl: 'https://api.openai.com/v1').credential()

        then:
        resolved == 'sk-real'
    }

    // -----------------------------------------------------------------------
    // A credential the endpoint gate WITHHELD is not a missing one. The placeholder exists for a
    // genuine no-credential local endpoint; substituting it for a misroute guarantees an opaque
    // 401 where a diagnosis was available (and, on pi, shadows an out-of-band key).
    // -----------------------------------------------------------------------

    def 'a withheld credential never becomes the placeholder'() {
        expect: 'the exact shape the gate produces: null apiKey, a baseUrl set, withheld flagged'
        new AgentRunnerRequest(model: 'openai/llama-3.3-70b', baseUrl: 'http://gw.corp/v1',
            credentialWithheld: true).credential() == null

        and: 'without the flag the very same request DOES get the placeholder -- that is the whole'
        // difference the flag carries, and why "nothing resolved" had to stop meaning both things
        new AgentRunnerRequest(model: 'openai/llama-3.3-70b', baseUrl: 'http://gw.corp/v1').credential() ==
            AgentRunnerRequest.PLACEHOLDER_API_KEY
    }

    def 'the withheld flag short-circuits D5 too, rather than reporting the wrong failure'() {
        when: 'a known provider host with a credential that resolved in ANOTHER namespace'
        // D5 would say "missing anthropic credential", which is false: one was found, for openai,
        // and refused. The runner that owns the decision reports the real cause.
        def result = new AgentRunnerRequest(model: 'openai/gpt-4o', baseUrl: 'https://api.anthropic.com/v1',
            credentialWithheld: true).credential()

        then:
        noExceptionThrown()
        result == null
    }

    def 'the withheld flag defaults to false, so an untouched request behaves exactly as before'() {
        expect:
        !new AgentRunnerRequest(model: 'openai/gpt-5-mini', apiKey: 'sk-real').credentialWithheld
        new AgentRunnerRequest(model: 'openai/gpt-5-mini', apiKey: 'sk-real').credential() == 'sk-real'
    }

    def 'the D5 error names agent.apiProvider, the third way out'() {
        when: 'the namespace this message reports was INFERRED from the host, so it is overridable'
        AgentRunnerRequest.credentialFor(null, 'https://api.anthropic.com/v1')

        then:
        def error = thrown(AbortOperationException)
        error.message.contains('`agent.apiKey`')
        error.message.contains('`NXF_AGENT_API_KEY`')
        error.message.contains('`agent.apiProvider`')
    }

    def 'the provider namespace travels on the request, unlike the credential'() {
        given: 'a runner reporting a missing credential must name the variables the ladder ACTUALLY'
        // consulted; that is a driver-side resolution, so it has to be carried
        def req = new AgentRunnerRequest(model: 'openai/gpt-4o', baseUrl: 'https://openrouter.ai/api/v1',
            apiProvider: 'openrouter', apiKey: 'sk-or')

        expect:
        req.apiProvider == 'openrouter'
        and: 'and it is not a secret, so unlike apiKey it is not excluded from toString()'
        req.toString().contains('openrouter')
        !req.toString().contains('sk-or')
    }
}
