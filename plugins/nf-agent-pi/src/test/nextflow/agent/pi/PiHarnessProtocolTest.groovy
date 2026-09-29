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

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.TimeUnit

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import nextflow.agent.AgentProtocolSpec
import nextflow.agent.AgentRunnerRequest
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Timeout

/**
 * Drives the REAL {@code harness/runner.mjs} over its JSONL protocol, in the role the
 * {@code agent-rpc} proxy plays inside the runner image.
 *
 * <p>The harness is the largest piece of the runtime that moved out of the plugin jar and into
 * the container image, and the image is not built by this project's Gradle build -- so without
 * this spec nothing in the build exercises it at all. The Pi SDK it imports is replaced by a
 * scripted stub installed into a throwaway {@code node_modules} (see
 * {@code src/testResources/harness/pi-coding-agent-stub.mjs}), which keeps the build free of
 * {@code npm}, of the network and of provider credentials while still running the harness's own
 * code: its framing, its tool brokering, its structured-output contract and the session it builds.
 *
 * <p>The stub ships no builtin tools, and that is faithful: a runner-native tool is the SDK's, not
 * the harness's, so what the harness controls -- and all these specs can therefore assert -- is the
 * {@code tools:} ALLOWLIST it hands to {@code createAgentSession}.
 *
 * <p>What a stub cannot catch is a breaking change in the real
 * {@code @earendil-works/pi-coding-agent}: it pins the API SURFACE the harness depends on (the
 * factory names, the session event shapes, the tool result contract), so a rename on OUR side
 * fails here, while a rename on THEIRS surfaces when the image is built against a new pinned
 * version. Running the SDK for real needs {@code npm} and provider credentials, which is exactly
 * what moving the runtime into the image took out of this build.
 *
 * <p>Skipped when {@code node} is not on PATH; the image build is what actually needs it.
 */
@Timeout(value = 60, unit = TimeUnit.SECONDS)
@Requires({ PiHarnessProtocolTest.nodeAvailable() })
class PiHarnessProtocolTest extends Specification {

    /** The invocation identity the proxy stamps on every frame in both directions. */
    private static final String INVOCATION = 'inv-1'

    @TempDir Path folder

    private Harness harness
    private Path workDir

    static boolean nodeAvailable() {
        try {
            return new ProcessBuilder('node', '--version')
                .redirectOutput(ProcessBuilder.Redirect.DISCARD)
                .redirectError(ProcessBuilder.Redirect.DISCARD)
                .start()
                .waitFor() == 0
        }
        catch( Exception e ) {
            return false
        }
    }

    def cleanup() {
        harness?.close()
    }

    def 'should announce the protocol version before anything else'() {
        when: 'the harness is started, as the proxy starts it'
        harness = start()

        then: 'the ready frame the proxy waits for (agent-rpc/main.go) comes first, unprompted'
        harness.next() == [type: 'ready', protocolVersion: 2]
    }

    def 'should reject a start frame it cannot act on'() {
        given:
        harness = start()
        harness.next()

        when:
        harness.send([type: 'start', protocolVersion: version, invocationId: INVOCATION, spec: specOf(spec)])
        final frame = harness.nextFrame()

        then: 'the error names the invocation, so the driver can attribute it'
        frame.type == 'error'
        frame.invocationId == INVOCATION
        frame.message.contains(reason)

        where: 'the current version is 2, so 1 is an OLD driver and 3 a newer one - both refused'
        version | spec                            || reason
        1       | [:]                             || 'Unsupported protocol version'
        3       | [:]                             || 'Unsupported protocol version'
        2       | [model: 'gpt-4o']               || 'Invalid model identifier'
        2       | [model: 'openai/no-such-model'] || 'Unknown Pi model'
    }

    def 'should report a malformed protocol line instead of dying on it'() {
        given:
        harness = start()
        harness.next()

        when: 'a line the host never should have written'
        harness.sendRaw('{ this is not json')
        final frame = harness.nextFrame()

        then:
        frame.type == 'error'
        frame.code == 'invalid_json'
    }

    def 'should complete with the assistant text, over a stdout that carries protocol frames only'() {
        given: 'a model that answers in one turn'
        harness = start([[[text: 'the capital is Paris']]])
        harness.next()

        when:
        harness.send(startFrame(specOf()))
        final frame = harness.nextFrame()

        // the stub logs through console.log while the run is in flight; the harness routes every
        // diagnostic channel to stderr, which is why each line read here parses as a frame
        then:
        frame == [
            type: 'complete',
            invocationId: INVOCATION,
            output: 'the capital is Paris',
            resolvedModel: 'openai/gpt-4o' ]
    }

    def 'should broker a declared tool to the host and feed the result back to the model'() {
        given: 'the model calls the tool, then answers with what the tool returned'
        harness = start([[[tool: 'word_stats', args: [path: 'reads.txt'], echo: true]]])
        harness.next()

        when:
        harness.send(startFrame(specOf(toolSpecs: [toolSpec('word_stats')])))
        final call = harness.nextFrame()

        then: 'a Nextflow tool is a process on the DRIVER - the harness asks the host to run it'
        call.type == 'tool_call'
        call.invocationId == INVOCATION
        call.name == 'word_stats'
        call.arguments == [path: 'reads.txt']
        and: 'and traces the dispatch for `-with-agent-trace`'
        harness.traces*.subMap(['event', 'name']) == [[event: 'tool_start', name: 'word_stats']]

        when: 'the host replies with the process output'
        harness.send([type: 'tool_result', invocationId: INVOCATION, callId: call.callId, result: '{"words":42}'])
        final done = harness.nextFrame()

        then: 'the model saw the result - it is echoed back as the answer'
        done.type == 'complete'
        done.output == '{"words":42}'
        and:
        harness.traces*.event == ['tool_start', 'tool_end']
    }

    def 'should refuse a tool result stamped with another invocation'() {
        given:
        harness = start([[[tool: 'word_stats', args: [:]]]])
        harness.next()
        harness.send(startFrame(specOf(toolSpecs: [toolSpec('word_stats')])))
        final call = harness.nextFrame()

        when: 'a result arrives for a different invocation than the one in flight'
        harness.send([type: 'tool_result', invocationId: 'someone-else', callId: call.callId, result: 'x'])
        final frame = harness.nextFrame()

        then:
        frame.type == 'error'
        frame.message.contains('Mismatched invocationId')
    }

    def 'should terminate on final_answer and return its arguments as the output'() {
        given: 'structured output is requested, and the model calls the schema-bound tool'
        harness = start([[[tool: 'final_answer', args: [capital: 'Paris', country: 'France']]]])
        harness.next()

        when:
        harness.send(startFrame(specOf(outputSchema: CAPITAL_SCHEMA)))
        final frame = harness.nextFrame()

        then:
        frame.type == 'complete'
        new JsonSlurper().parseText(frame.output as String) == [capital: 'Paris', country: 'France']
    }

    def 'should take one corrective turn when the model answers structured output with prose'() {
        given: 'turn one is ordinary text; turn two is the re-prompt the harness issues'
        harness = start([
            [[text: 'The capital of France is Paris.']],
            [[tool: 'final_answer', args: [capital: 'Paris', country: 'France']]] ])
        harness.next()

        when:
        harness.send(startFrame(specOf(outputSchema: CAPITAL_SCHEMA)))
        final frame = harness.nextFrame()

        then: 'the prose is discarded - only the schema-bound arguments are the output'
        frame.type == 'complete'
        new JsonSlurper().parseText(frame.output as String) == [capital: 'Paris', country: 'France']
    }

    def 'should attribute an empty answer to the provider failure that caused it'() {
        given: 'the SDK reports a failed model call as an assistant message with no content'
        harness = start([[[providerError: 'HTTP 500 (request id req_abc123)']]])
        harness.next()

        when:
        harness.send(startFrame(specOf()))
        final frame = harness.nextFrame()

        then: 'an exhausted retry chain is not misreported as the model declining to answer'
        frame.type == 'error'
        frame.message.contains('Pi returned no final assistant text')
        frame.message.contains('req_abc123')
    }

    // --- the runner split: brokered descriptors vs runner-native names

    def 'should broker EVERY declared descriptor, whatever it is named'() {
        given: 'a tool named exactly like the capability the harness used to serve in-process'
        harness = start([[[tool: 'filesystem', args: [path: 'reads.txt'], echo: true]]])
        harness.next()

        when:
        harness.send(startFrame(specOf(toolSpecs: [toolSpec('filesystem')])))
        final call = harness.nextFrame()

        then: 'there is no local branch left to fall into - a descriptor means the DRIVER owns it'
        call.type == 'tool_call'
        call.name == 'filesystem'

        when:
        harness.send([type: 'tool_result', invocationId: INVOCATION, callId: call.callId, result: 'from the driver'])
        final done = harness.nextFrame()

        then:
        done.type == 'complete'
        done.output == 'from the driver'
    }

    def 'should enable the runner-native tools through the session allowlist, beside the brokered ones'() {
        given: 'the answer is the session the harness built'
        harness = start([[[runtimeState: true]]])
        harness.next()

        when: 'the agent selected one `nf:module_run` tool and four `fs:`/`shell:` leaves'
        harness.send(startFrame(specOf(
            toolSpecs: [toolSpec('word_stats')],
            nativeToolNames: ['read', 'write', 'grep', 'bash'] )))
        final session = new JsonSlurper().parseText(harness.nextFrame().output as String).session

        then: 'a native name enters the ALLOWLIST, which is how an SDK builtin is turned on;'
        // it is never a customTool, because the runner - not this harness - implements it
        session.tools == ['word_stats', 'read', 'write', 'grep', 'bash']
        session.customTools == ['word_stats']

        and: 'the allowlist is the whole gate, so `noTools` is gone rather than contradicting it'
        session.noTools == null

        and: 'the builtins are rooted at the task work dir - that is the pi-side sandbox'
        session.cwd == workDir.toString()
    }

    def 'should enable nothing at all when the agent declared no tools'() {
        given:
        harness = start([[[runtimeState: true]]])
        harness.next()

        when: 'G7: no `tools` directive, hence no descriptors and no native names'
        harness.send(startFrame(specOf()))
        final session = new JsonSlurper().parseText(harness.nextFrame().output as String).session

        then: 'an EMPTY allowlist, not an absent one - pi would otherwise enable its four defaults'
        session.tools == []
    }

    def 'should stop a runaway agent at maxIterations'() {
        given: 'two tool turns against a budget of one'
        harness = start([[
            [tool: 'word_stats', args: [:]],
            [tool: 'word_stats', args: [:]] ]])
        harness.next()

        when:
        harness.send(startFrame(specOf(maxIterations: 1, toolSpecs: [toolSpec('word_stats')])))
        final call = harness.nextFrame()
        harness.send([type: 'tool_result', invocationId: INVOCATION, callId: call.callId, result: 'ok'])
        final frame = harness.nextFrame()

        then: 'the budget is counted on this side, so a second turn never reaches the host'
        frame.type == 'error'
        frame.message.contains('exceeded the maximum number of tool-call iterations (1)')
    }

    def 'should retarget the provider catalog at the endpoint resolved by the driver'() {
        given: 'the answer is whatever the harness asked of the ModelRuntime'
        harness = start([[[runtimeState: true]]])
        harness.next()

        when: '`agent.baseUrl` resolved to a compatible endpoint and travelled inside the spec'
        harness.send(startFrame(specOf(baseUrl: 'http://localhost:8000/v1')))
        final frame = harness.nextFrame()

        then: 'registerProvider rewrites the endpoint of the catalog, which is the only seam the'
        // SDK sanctions -- the endpoint must NOT be spread onto the resolved model, because the
        // catalog re-resolves Model instances by id and would drop it
        new JsonSlurper().parseText(frame.output as String).providers == [openai: [baseUrl: 'http://localhost:8000/v1']]
    }

    def 'should install no credential from the frame AgentRpcBroker actually builds'() {
        given:
        harness = start([[[runtimeState: true]]])
        harness.next()

        when: 'an endpoint, and deliberately no credential: a runtime key OWNS its provider, so'
        // one pushed from the driver would shadow exactly the credential the `env` scope or
        // `secret` directive delivered to this container. Pi is left to resolve its own.
        harness.send(startFrame(specOf(baseUrl: 'http://localhost:8000/v1')))
        final frame = harness.nextFrame()

        then:
        new JsonSlurper().parseText(frame.output as String).apiKeys == [:]
    }

    def 'should install a credential the start frame carries beside the spec'() {
        given: 'the receiving half, kept live for a sender that can establish the credential is'
        // for THIS provider -- see the send-site comment in AgentRpcBroker
        harness = start([[[runtimeState: true]]])
        harness.next()

        when: 'the credential travels BESIDE the portable spec, never inside it'
        harness.send(startFrame(specOf()) + [apiKey: 'gateway-credential-4a71b2'])
        final frame = harness.nextFrame()

        then: 'it is installed in memory for this process, never in its environment'
        new JsonSlurper().parseText(frame.output as String).apiKeys == [openai: 'gateway-credential-4a71b2']

        and: 'and the payload that crosses the wire never carried it'
        !specOf().containsKey('apiKey')
    }

    def 'should report an unexpected EOF when the host closes the protocol input'() {
        given:
        harness = start()
        harness.next()

        when: 'the proxy dies, or the driver goes away, before the run finishes'
        harness.closeInput()
        final frame = harness.nextFrame()

        then: 'the container does not exit silently with no answer and no reason'
        frame.type == 'error'
        frame.code == 'unexpected_eof'
    }

    // --- fixtures

    private static final Map CAPITAL_SCHEMA = [
        type: 'object',
        properties: [capital: [type: 'string'], country: [type: 'string']],
        required: ['capital', 'country'],
        additionalProperties: false ]

    private static Map toolSpec(String name) {
        return [
            name: name,
            description: "Nextflow tool ${name}".toString(),
            inputSchema: [type: 'object', properties: [:], additionalProperties: true] ]
    }

    /**
     * The request payload the broker sends, built by the SAME core class the broker uses
     * ({@link AgentProtocolSpec}), so a spec field renamed there fails here rather than in the
     * container.
     */
    private Map specOf(Map overrides = [:]) {
        final defaults = [
            model: 'openai/gpt-4o',
            prompt: 'What is the capital of France?',
            maxIterations: 10,
            workDir: workDir.toString() ]
        return AgentProtocolSpec.fromRequest(new AgentRunnerRequest(defaults + overrides))
    }

    private static Map startFrame(Map spec) {
        return [type: 'start', protocolVersion: 2, invocationId: INVOCATION, spec: spec]
    }

    /**
     * Start the real harness against the stub SDK. Both are COPIED into a throwaway directory:
     * Node resolves a bare import by walking up from the importing file, so running the harness
     * where it lives would pick up a developer's `npm ci` tree instead of the stub -- and the
     * point of this build is that no such tree has to exist.
     */
    private Harness start(List plan = []) {
        final sandbox = Files.createDirectories(folder.resolve('sandbox'))
        workDir = Files.createDirectories(folder.resolve('work'))
        Files.copy(harnessSource(), sandbox.resolve('runner.mjs'))
        final module = Files.createDirectories(sandbox.resolve('node_modules/@earendil-works/pi-coding-agent'))
        module.resolve('package.json').text = JsonOutput.toJson([
            name: '@earendil-works/pi-coding-agent',
            version: '0.0.0-stub',
            type: 'module',
            exports: './index.mjs' ])
        module.resolve('index.mjs').text = stubSource()
        return new Harness(sandbox, JsonOutput.toJson(plan))
    }

    /** The harness under test, handed over by the Gradle test task -- see build.gradle. */
    private static Path harnessSource() {
        final path = System.getProperty('nf.agent.pi.harness')
        assert path, 'nf.agent.pi.harness is unset - run this spec through Gradle (:plugins:nf-agent-pi:test)'
        final file = Path.of(path)
        assert Files.exists(file), "the harness was not found: ${file}"
        return file
    }

    private static String stubSource() {
        final url = PiHarnessProtocolTest.getResource('/harness/pi-coding-agent-stub.mjs')
        assert url, 'the Pi SDK stub is missing from the test resources'
        return url.text
    }

    /**
     * The proxy's half of the protocol: JSON frames in on stdin, JSON frames out on stdout, one
     * per line. Mirrors what agent-rpc/main.go does between the driver's gRPC stream and the
     * harness process.
     */
    private static class Harness implements Closeable {

        private final Process process
        private final BufferedReader output
        private final Writer input

        /** Trace frames seen so far; they interleave with the frames a spec asserts on. */
        final List<Map> traces = []

        Harness(Path dir, String plan) {
            final builder = new ProcessBuilder('node', 'runner.mjs')
                .directory(dir.toFile())
                .redirectError(ProcessBuilder.Redirect.INHERIT)
            builder.environment().put('PI_STUB_PLAN', plan)
            process = builder.start()
            output = new BufferedReader(new InputStreamReader(process.getInputStream(), 'UTF-8'))
            input = new OutputStreamWriter(process.getOutputStream(), 'UTF-8')
        }

        void send(Map frame) {
            sendRaw(JsonOutput.toJson(frame))
        }

        void sendRaw(String line) {
            input.write(line)
            input.write('\n')
            input.flush()
        }

        void closeInput() {
            input.close()
        }

        /** The next frame, or {@code null} at end of stream. */
        Map next() {
            final line = output.readLine()
            return line != null ? new JsonSlurper().parseText(line) as Map : null
        }

        /** The next non-trace frame; traces are collected into {@link #traces} on the way. */
        Map nextFrame() {
            Map frame = next()
            while( frame?.type == 'trace' ) {
                traces.add(frame)
                frame = next()
            }
            return frame
        }

        @Override
        void close() {
            try { input.close() } catch( IOException e ) { /* already closed */ }
            process.destroy()
            process.waitFor(10, TimeUnit.SECONDS)
        }
    }
}
