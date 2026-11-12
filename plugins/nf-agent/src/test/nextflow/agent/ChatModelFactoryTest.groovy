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

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.json.JsonSchema
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.Timeout

class ChatModelFactoryTest extends Specification {

    def 'should split provider and model id'() {
        expect:
        ChatModelFactory.providerOf('openai/gpt-5-mini') == 'openai'
        ChatModelFactory.modelOf('openai/gpt-5-mini') == 'gpt-5-mini'
    }

    def 'the prefix this factory reads is the same one core resolves credentials from'() {
        given: 'core lower-cases the prefix and this factory did not, so a mixed-case id made them'
        // disagree about whose key was in flight -- core resolved and TRANSMITTED an OPENAI_API_KEY
        // for a model this factory then rejected as an unsupported provider
        expect:
        ChatModelFactory.providerOf(MODEL) == AgentConfig.providerPrefixOf(MODEL)

        where:
        MODEL << ['openai/gpt-5-mini', 'OpenAI/gpt-4o', 'OPENAI/gpt-4o', 'Anthropic/claude-sonnet-4']
    }

    def 'a mixed-case openai prefix builds a model instead of being refused'() {
        when: 'the direct consequence of the agreement above'
        ChatModel model = new ChatModelFactory().createModel('OpenAI/gpt-4o', 30, null, (Double) null, 'sk-test', null, null, false)
        then:
        model != null
    }

    def 'should fail for a model id without a provider'() {
        when:
        ChatModelFactory.providerOf('gpt-5-mini')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should fail for an unknown provider'() {
        when:
        new ChatModelFactory().createModel('acme/whatever', 30, null, (Double) null, 'sk-test', null, null, false)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.toLowerCase().contains('provider')
        and: 'the message says it is the WIRE PROTOCOL that is unsupported, and points at baseUrl'
        e.message.contains('wire protocol')
        e.message.contains('agent.baseUrl')
    }

    def 'should fail when neither a credential nor an endpoint resolves'() {
        when:
        new ChatModelFactory().createModel('openai/gpt-5-mini', 30, null, (Double) null, null, null, null, false)
        then:
        def e = thrown(IllegalArgumentException)
        and: 'the message names the whole resolution ladder, not one variable'
        e.message.contains('agent.apiKey')
        e.message.contains('NXF_AGENT_API_KEY')
        e.message.contains('OPENAI_API_KEY')
        and: 'plus the one rung this factory cannot see -- the namespace an explicit option redirected to'
        // credentialFor has already aborted for an endpoint whose HOST names a provider, so the
        // only way the ladder read something other than the OPENAI_ variables and still arrived
        // here empty is an `agent.apiProvider` the request does not carry
        e.message.contains('agent.apiProvider')
    }

    def 'should refuse to present a placeholder to a well-known provider endpoint'() {
        given: 'design D5. The placeholder assumes the endpoint needs no credential -- true of a'
        // local vLLM/Ollama, plainly false of a provider's own API. Substituting it there buys an
        // opaque 401 one request later instead of a diagnosis now.
        when:
        new ChatModelFactory().createModel('openai/gpt-4o', 30, null, (Double) null, null, ENDPOINT, null, false)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains(PROVIDER)
        e.message.contains(ENDPOINT)
        e.message.contains(VAR)

        where:
        ENDPOINT                       | PROVIDER     | VAR
        'https://api.openai.com/v1'    | 'openai'     | 'OPENAI_API_KEY'
        'https://openrouter.ai/api/v1' | 'openrouter' | 'OPENROUTER_API_KEY'
    }

    @Timeout(30)
    def 'the widened credential ladder reaches the builder for a non-openai namespace'() {
        given: 'the whole point of D1/D2: `openai` is the WIRE PROTOCOL, and the credential comes'
        // from whichever namespace `agent.apiProvider` names -- OpenRouter here, which the openai
        // carve-out had no spelling for at all. Core resolves it; this factory only presents it.
        def authorizations = Collections.synchronizedList([])
        def server = HttpServer.create(new InetSocketAddress(0), 0)
        server.createContext('/', new HttpHandler() {
            @Override
            void handle(HttpExchange exchange) {
                authorizations.add(exchange.requestHeaders.getFirst('Authorization'))
                final body = '{"id":"1","object":"chat.completion","created":0,"model":"local","choices":[{"index":0,"message":{"role":"assistant","content":"pong"},"finish_reason":"stop"}]}'
                exchange.responseHeaders.add('Content-Type', 'application/json')
                exchange.sendResponseHeaders(200, body.bytes.length)
                exchange.responseBody.withCloseable { it.write(body.bytes) }
            }
        })
        server.start()

        and: 'the endpoint stands in for the OpenRouter gateway, vouched for by an explicit namespace'
        def endpoint = "http://localhost:${server.address.port}/v1".toString()
        def config = new AgentConfig([apiProvider: 'openrouter', baseUrl: endpoint],
            [OPENROUTER_API_KEY: 'sk-or-resolved', OPENAI_API_KEY: 'sk-openai-must-not-travel'])

        when: 'the pair core resolved is handed to the factory, as LangChainAgentRunner does'
        def model = new ChatModelFactory().createModel('openai/gpt-4o', 30, null, (Double) null,
            config.apiKeyFor('openai/gpt-4o'), config.baseUrlFor('openai/gpt-4o'),
            config.apiProviderFor('openai/gpt-4o'), config.credentialWithheldFor('openai/gpt-4o'))
        def answer = model.chat('ping')

        then: 'the OpenRouter credential is what was presented, and the OpenAI one never left the driver'
        answer == 'pong'
        authorizations == ['Bearer sk-or-resolved']

        cleanup:
        server?.stop(0)
    }

    def 'should build a model with an endpoint and no credential (no network call)'() {
        when: 'a local endpoint needing no credential - a placeholder is sent instead of failing'
        ChatModel model = new ChatModelFactory().createModel('openai/llama-3.3-70b', 30, null, (Double) null, null, 'http://localhost:8000/v1', null, false)
        then:
        noExceptionThrown()
        model != null
    }

    def 'should build an openai model when api key is present (no network call)'() {
        given:
        def factory = new ChatModelFactory()
        when:
        ChatModel model = factory.createModel('openai/gpt-5-mini', 30, null, (Double) null, 'sk-test', null, null, false)
        then:
        model != null
    }

    @Timeout(30)
    def 'should send the chat request to the configured endpoint'() {
        given: 'a loopback stub standing in for an OpenAI-compatible endpoint (vLLM, Ollama, a gateway)'
        String requestedPath = null
        String authorization = null
        def server = HttpServer.create(new InetSocketAddress(0), 0)
        server.createContext('/', new HttpHandler() {
            @Override
            void handle(HttpExchange exchange) {
                requestedPath = exchange.requestURI.path
                authorization = exchange.requestHeaders.getFirst('Authorization')
                final body = '{"id":"1","object":"chat.completion","created":0,"model":"local","choices":[{"index":0,"message":{"role":"assistant","content":"pong"},"finish_reason":"stop"}]}'
                exchange.responseHeaders.add('Content-Type', 'application/json')
                exchange.sendResponseHeaders(200, body.bytes.length)
                exchange.responseBody.withCloseable { it.write(body.bytes) }
            }
        })
        server.start()

        and: 'no credential resolved, so the placeholder is sent (D8)'
        def model = new ChatModelFactory().createModel('openai/llama-3.3-70b', 30, null, (Double) null, null, "http://localhost:${server.address.port}/v1".toString(), null, false)

        when:
        def answer = model.chat('ping')

        then: 'the request went to the configured endpoint, not to api.openai.com'
        answer == 'pong'
        requestedPath == '/v1/chat/completions'
        and: 'carrying the placeholder SHARED with the pi runner, so the two cannot drift'
        authorization == "Bearer ${AgentRunnerRequest.PLACEHOLDER_API_KEY}"

        cleanup:
        server?.stop(0)
    }

    @Timeout(30)
    def 'one shared factory serves two different credentials (per-call, not per-instance)'() {
        given: 'LangChainAgentRunner holds ONE long-lived factory shared by concurrently executing'
        // agent tasks, so a per-request credential kept as a FIELD would be a data race across
        // parallel agents -- and a field initialized from the environment would put a second,
        // divergent resolution ladder next to the core one (design D4)
        def authorizations = Collections.synchronizedList([])
        def server = HttpServer.create(new InetSocketAddress(0), 0)
        server.createContext('/', new HttpHandler() {
            @Override
            void handle(HttpExchange exchange) {
                authorizations.add(exchange.requestHeaders.getFirst('Authorization'))
                final body = '{"id":"1","object":"chat.completion","created":0,"model":"local","choices":[{"index":0,"message":{"role":"assistant","content":"pong"},"finish_reason":"stop"}]}'
                exchange.responseHeaders.add('Content-Type', 'application/json')
                exchange.sendResponseHeaders(200, body.bytes.length)
                exchange.responseBody.withCloseable { it.write(body.bytes) }
            }
        })
        server.start()
        and:
        def endpoint = "http://localhost:${server.address.port}/v1".toString()
        def factory = new ChatModelFactory()

        when: 'the SAME factory builds two models with different credentials'
        factory.createModel('openai/gpt-5-mini', 30, null, (Double) null, 'sk-alpha', endpoint, null, false).chat('ping')
        factory.createModel('openai/gpt-5-mini', 30, null, (Double) null, 'sk-beta', endpoint, null, false).chat('ping')

        then: 'each request carried its own key -- no cross-talk, and no placeholder substituted'
        authorizations == ['Bearer sk-alpha', 'Bearer sk-beta']

        cleanup:
        server?.stop(0)
    }

    def 'should build an openai model with an explicit temperature (no network call)'() {
        given:
        def factory = new ChatModelFactory()
        when: 'a pinned temperature is applied on the builder (0.0 must not throw)'
        ChatModel model = factory.createModel('openai/gpt-4o', 30, null, 0.0d, 'sk-test', null, null, false)
        then:
        model != null
    }

    def 'a WITHHELD credential is a distinct, fatal, named failure -- never a placeholder'() {
        given: 'the driver resolved OPENAI_API_KEY and the endpoint gate refused to send it to a'
        // gateway OpenAI does not own. Substituting the D8 placeholder there guarantees an opaque
        // 401 where a diagnosis is available; langchain4j has no other credential source, so it is
        // fatal here -- deliberately UNLIKE the pi runner, which only warns.
        def config = new AgentConfig([baseUrl: 'https://gw.corp/v1'], [OPENAI_API_KEY: 'sk-oai'])

        when:
        new ChatModelFactory().createModel('openai/gpt-4o', 30, null, (Double) null,
            config.apiKeyFor('openai/gpt-4o'), config.baseUrlFor('openai/gpt-4o'),
            config.apiProviderFor('openai/gpt-4o'), config.credentialWithheldFor('openai/gpt-4o'))

        then:
        def e = thrown(IllegalArgumentException)
        and: 'it names the three ways out, and the endpoint the credential was refused for'
        e.message.contains('`agent.apiKey`')
        e.message.contains('NXF_AGENT_API_KEY')
        e.message.contains('agent.apiProvider')
        e.message.contains('https://gw.corp/v1')
        and: 'and it is NOT the generic "missing credential" message: one was found'
        !e.message.contains('Missing LLM provider credential')
    }

    def 'the withheld error names the default endpoint when no baseUrl resolved'() {
        given: 'the F1 case: an apiProvider that is not the model prefix, so the request would go'
        // to the prefix provider's own default endpoint with somebody else's key
        def config = new AgentConfig([apiProvider: 'openrouter'], [OPENROUTER_API_KEY: 'sk-or'])

        when:
        new ChatModelFactory().createModel('openai/gpt-4o', 30, null, (Double) null,
            config.apiKeyFor('openai/gpt-4o'), config.baseUrlFor('openai/gpt-4o'),
            config.apiProviderFor('openai/gpt-4o'), config.credentialWithheldFor('openai/gpt-4o'))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('`openrouter`')
        e.message.contains('default `openai` endpoint')
    }

    def 'the missing-credential message names the RESOLVED namespace, not always OpenAI'() {
        when: 'agent.apiProvider redirected the namespace and nothing was exported for it'
        // the old message asserted OPENAI_API_KEY unconditionally -- the very thing its own comment
        // claimed to avoid
        new ChatModelFactory().createModel('openai/gpt-4o', 30, null, (Double) null, null, null, 'mistral', false)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('MISTRAL_API_KEY')
        e.message.contains('NXF_AGENT_API_KEY')
        e.message.contains('agent.apiProvider')
        and: 'the variable that was NOT read is not suggested'
        !e.message.contains('OPENAI_API_KEY')
    }

    def 'should build an openai model with a structured-output schema (no network call)'() {
        given:
        def factory = new ChatModelFactory()
        def schema = JsonSchemaMapper.toJsonSchema('Answer', [
            type: 'object',
            properties: [answer: [type: 'string']],
            required: ['answer'],
            additionalProperties: false,
        ])
        when:
        ChatModel model = factory.createModel('openai/gpt-5-mini', 30, schema, (Double) null, 'sk-test', null, null, false)
        then:
        model != null
    }
}
