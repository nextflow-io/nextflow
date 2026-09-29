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
package nextflow.script

import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.spi.FileSystemProvider

import groovy.json.JsonOutput
import nextflow.SysEnv
import nextflow.agent.AgentOutputMode
import nextflow.agent.AgentOutputPlan
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.SkillDescriptor
import nextflow.agent.SkillResource
import nextflow.agent.ToolDescriptor
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHolder
import nextflow.processor.TaskPath
import nextflow.script.TokenValRef
import spock.lang.Specification
import nextflow.agent.rpc.AgentRpcConfig

class AgentDefTest extends Specification {

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

    private AgentDef makeAgent(BaseScript script, String name) {
        def prompt = new PromptDef({ -> 'hello' }, 'hello')
        return new AgentDef(script, name, [:], [], [], prompt)
    }

    def 'should construct an AgentDef with name and content'() {
        given:
        def script = Mock(BaseScript)

        when:
        def agent = makeAgent(script, 'eval_agent')

        then:
        agent.name == 'eval_agent'
        agent.simpleName == 'eval_agent'
        agent.type == 'agent'
    }

    def 'should preserve a remote path scheme in agent input json'() {
        given:
        final provider = Mock(FileSystemProvider)
        final fs = Mock(FileSystem)
        final path = Mock(Path)
        path.getFileSystem() >> fs
        fs.provider() >> provider
        provider.getScheme() >> 'mock'
        path.toAbsolutePath() >> path
        path.toUri() >> URI.create('mock://bucket/reads.fastq')

        expect:
        AgentDef.toJson([reads: path]) == '{"reads":"mock://bucket/reads.fastq"}'
        and: 'nesting it in a Map-typed input -- a shape compile-time inference never stages -- keeps the scheme too'
        AgentDef.toJson([refs: [primary: path]]) == '{"refs":{"primary":"mock://bucket/reads.fastq"}}'
    }

    def 'should render a STAGED input by its work-dir stage name'() {
        given: 'the value shape staging produces: a work-dir-relative view of an input file'
        final source = Path.of('/data/run1/contigs.fa')
        final staged = new TaskPath(new FileHolder(source))

        expect: 'the JSON carries the name the model can open inside the task dir'
        AgentDef.toJson(staged) == '"contigs.fa"'
        and: 'which is exactly what the prompt interpolation of the same input produces --'
        // this is the invariant: one input, one rendering, in the prompt and in the JSON
        AgentDef.toJson(staged) == JsonOutput.toJson("${staged}".toString())
        and: 'a staged path nested inside a record input renders the same way'
        AgentDef.toJson([id: 's1', seq: staged]) == '{"id":"s1","seq":"contigs.fa"}'
    }

    def 'should fail on run() when a tool-free agent declares zero outputs (task path)'() {
        given:
        // no tools/skills -> task path, where the one-input guard is gone; a
        // zero-output agent is still an error
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo') // no inputs, no outputs, no tools

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare exactly one output')
    }

    def 'should fail on run() with an input arity mismatch (task path)'() {
        given:
        // 1 declared input, 1 declared output, no tools -> task path; called with 0 args
        def script = Mock(BaseScript)
        def inp = new AgentBuilder.AgentInput('q', String)
        def out = new AgentBuilder.AgentOutput('a', String)
        def agent = new AgentDef(script, 'foo', [:] as Map<String,Object>, [inp], [out], new PromptDef({ -> 'h' }, 'h'))

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('expects 1 input channel')
    }

    def 'a tools agent lowers to the task path (legacy one-input guard removed)'() {
        given:
        // M-Tools: declaring `tools` no longer forces the legacy operator path — every agent
        // lowers to the task path, so the legacy exactly-one-input guard is gone (multiple
        // inputs are now permitted). The generalized guards apply instead; here the zero-output
        // guard fires (0 inputs / 0 args passes the arity check), proving task-path routing.
        def script = Mock(BaseScript)
        def agent = new AgentDef(script, 'foo', [tools: 'nf:module_run:someProc'] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h'))

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare exactly one output')
    }

    def 'should clone with a new name'() {
        given:
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo')

        when:
        def renamed = agent.cloneWithName('bar')

        then:
        renamed instanceof AgentDef
        renamed.name == 'bar'
        agent.name == 'foo' // original untouched
    }

    def 'cloneWithName preserves the declared base name'() {
        given:
        def agent = makeAgent(Mock(BaseScript), 'critic')

        when: 'the agent is included under an alias, then invoked in a named workflow scope'
        def aliased = (AgentDef) agent.cloneWithName('reviewer')
        def qualified = (AgentDef) agent.cloneWithName('WF:critic')

        then: 'the declared name survives, so `withName:critic` keeps matching'
        aliased.baseName == 'critic'
        aliased.name == 'reviewer'
        aliased.simpleName == 'reviewer'
        and:
        qualified.baseName == 'critic'
        qualified.name == 'WF:critic'
        qualified.simpleName == 'critic'
    }

    def 'should expose the declared labels'() {
        expect:
        new AgentDef(Mock(BaseScript), 'a', [label: ['big', 'fast']] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).labels == ['big', 'fast']
        and: 'an agent with no label declaration has none'
        makeAgent(Mock(BaseScript), 'a').labels == []
    }

    def 'should expose the goal directive'() {
        given:
        def directives = [model: 'openai/gpt-5-mini', instruction: 'be careful', goal: 'assemble then QC'] as Map<String,Object>
        def prompt = new PromptDef({ -> 'hi' }, 'hi')
        def agent = new AgentDef(Mock(BaseScript), 'a', directives, [], [], prompt)

        expect:
        agent.goal == 'assemble then QC'
        agent.instruction == 'be careful'
    }

    def 'should return null goal when not declared'() {
        given:
        def agent = new AgentDef(Mock(BaseScript), 'a', [model: 'openai/gpt-5-mini'] as Map<String,Object>, [], [],
            new PromptDef({ -> 'hi' }, 'hi'))
        expect:
        agent.goal == null
    }

    static class TestRec implements nextflow.script.types.Record {}

    static class WrapPlan implements nextflow.script.types.Record {
        String title
        Long count
    }

    def 'buildWrapperSchema builds an object-root wrapper with per-output fragments'() {
        given:
        def outs = [
            new AgentBuilder.AgentOutput('rec', WrapPlan),
            new AgentBuilder.AgentOutput('n', Long),
            new AgentBuilder.AgentOutput('score', Double),
            new AgentBuilder.AgentOutput('flag', Boolean),
            new AgentBuilder.AgentOutput('label', String),
        ]

        when:
        def schema = AgentDef.buildWrapperSchema('agentX', outs)

        then:
        schema.type == 'object'
        schema.additionalProperties == false
        schema.required == ['rec', 'n', 'score', 'flag', 'label']
        (schema.properties.keySet() as List) == ['rec', 'n', 'score', 'flag', 'label']

        and: 'scalar fragments map to the right JSON-schema type'
        schema.properties.n.type == 'integer'
        schema.properties.score.type == 'number'
        schema.properties.flag.type == 'boolean'
        schema.properties.label.type == 'string'

        and: 'nested record fragment recursion is intact'
        schema.properties.rec.type == 'object'
        schema.properties.rec.properties.title.type == 'string'
        schema.properties.rec.properties.count.type == 'integer'
    }

    def 'buildWrapperSchema rejects an unsupported top-level output type'() {
        when:
        AgentDef.buildWrapperSchema('agentX', [new AgentBuilder.AgentOutput('p', java.nio.file.Path)])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('unsupported type')
        e.message.contains('`p`')
        and: 'the message enumerates the supported output set (plan §4.5/§9)'
        e.message.contains('supported:')
        e.message.contains('record type')
    }

    def 'scalarOutputSchema represents a Path as an exact string field'() {
        when:
        def schema = AgentDef.scalarOutputSchema(new AgentBuilder.AgentOutput('assembly_path', java.nio.file.Path))

        then:
        schema.type == 'object'
        schema.properties.assembly_path.type == 'string'
        schema.required == ['assembly_path']
        schema.additionalProperties == false
    }

    def 'decodeCanonicalOutput unwraps a scalar final answer contract'() {
        given:
        def stdout = '{"type":"complete","output":"{\\"answer\\":\\"HELLO\\"}"}'

        expect:
        new AgentOutputPlan(AgentOutputMode.SCALAR_CONTRACT, null).decode(stdout, 'answer', String) == 'HELLO'
    }

    // -----------------------------------------------------------------------
    // M2 resume: canonical BodyDef.source, key-sorted schema JSON, floating alias,
    // PromptDef valRefs (design D2/D3/D5)
    // -----------------------------------------------------------------------

    private AgentDef agentWith(Map directives, PromptDef prompt) {
        return new AgentDef(Mock(BaseScript), 'a', directives as Map<String,Object>, [], [], prompt)
    }

    def 'canonicalAgentSource is deterministic for the same effective inputs'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']], required: ['x'], additionalProperties: false]
        def a = agentWith([instruction: 'be careful', goal: 'do it'], new PromptDef({ -> 'p' }, 'the-prompt'))
        def b = agentWith([instruction: 'be careful', goal: 'do it'], new PromptDef({ -> 'p' }, 'the-prompt'))

        expect:
        a.canonicalAgentSource('openai/gpt-4o', 20, schema) == b.canonicalAgentSource('openai/gpt-4o', 20, schema)
    }

    def 'canonicalAgentSource includes the effective model (differs when model differs)'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']]]
        def a = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))

        expect: 'the effective model id (resolved this.model ?: default) is folded in'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema) != a.canonicalAgentSource('openai/gpt-4o-mini', 20, schema)

        and: 'the temperature line is the stable literal (task path leaves it unset)'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema).contains('temperature=default')
        !a.canonicalAgentSource('openai/gpt-4o', 20, schema).contains('temperature=0')
    }

    def 'canonicalAgentSource includes an explicitly selected runner'() {
        given:
        def agent = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))

        expect:
        agent.canonicalAgentSource('openai/gpt-4o', 20, null, null, 'pi') !=
            agent.canonicalAgentSource('openai/gpt-4o', 20, null, null, 'langchain4j')
        agent.canonicalAgentSource('openai/gpt-4o', 20, null, null, 'pi').startsWith('agentRunner=pi\n')
    }

    // -----------------------------------------------------------------------
    // M-Skills: skill identity folded into the resume cache key (the 4-arg
    // canonicalAgentSource overload). These are RED until the overload exists:
    // the reflective lookup below throws NoSuchMethodException on current source.
    // -----------------------------------------------------------------------

    /**
     * Invoke the (yet-to-exist) 4-arg {@code canonicalAgentSource(String, int, Map, List)} overload
     * reflectively so the test compiles against current source (no such method) and fails at RUNTIME
     * (RED) rather than breaking compilation of the whole test module — mirroring the reflective
     * private-method helpers used elsewhere in this class. Once the overload lands this resolves and
     * the assertions run (GREEN).
     */
    private static String canonical4(AgentDef agent, String model, int maxIter, Map schema, List<SkillDescriptor> skills) {
        def m = AgentDef.getDeclaredMethod('canonicalAgentSource', String, Integer.TYPE, Map, List)
        m.accessible = true
        return (String) m.invoke(agent, model, maxIter, schema, skills)
    }

    private static SkillDescriptor skill(String name, String content, List<SkillResource> resources = []) {
        return new SkillDescriptor(name, "desc of ${name}".toString(), content, resources)
    }

    def 'canonicalAgentSource 4-arg with null/empty skills equals the 3-arg form byte-for-byte'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']], required: ['x'], additionalProperties: false]
        def a = agentWith([instruction: 'be careful', goal: 'do it'], new PromptDef({ -> 'p' }, 'the-prompt'))
        def three = a.canonicalAgentSource('openai/gpt-4o', 20, schema)

        expect: 'a tool-free/skill-free agent keeps an IDENTICAL cache key (no spurious resume invalidation)'
        canonical4(a, 'openai/gpt-4o', 20, schema, null) == three
        canonical4(a, 'openai/gpt-4o', 20, schema, [] as List<SkillDescriptor>) == three
    }

    def 'canonicalAgentSource 4-arg differs once skills are present and when a skill changes'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']]]
        def a = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))
        def three = a.canonicalAgentSource('openai/gpt-4o', 20, schema)
        def withGreet = canonical4(a, 'openai/gpt-4o', 20, schema, [ skill('greet', 'v1 instructions') ])

        expect: 'declaring a skill changes the key (skills are a captured input now)'
        withGreet != three

        and: 'a changed SKILL.md body invalidates the cache (different content -> different key)'
        canonical4(a, 'openai/gpt-4o', 20, schema, [ skill('greet', 'v2 instructions') ]) != withGreet

        and: 'a changed bundled resource also invalidates the cache'
        canonical4(a, 'openai/gpt-4o', 20, schema, [ skill('greet', 'v1 instructions', [ new SkillResource('references/a.txt', 'AAA') ]) ]) != withGreet

        and: 'the fingerprint is order-independent (sorted by name) - same set of skills, reordered, same key'
        def s1 = skill('alpha', 'a-body')
        def s2 = skill('beta', 'b-body')
        canonical4(a, 'openai/gpt-4o', 20, schema, [ s1, s2 ]) == canonical4(a, 'openai/gpt-4o', 20, schema, [ s2, s1 ])
    }

    // -----------------------------------------------------------------------
    // M-Tools resume: tool identity folded into the resume cache key, which is what lets a
    // module/process-tool agent be cacheable at all (instead of a blanket `cache false`).
    // -----------------------------------------------------------------------

    private static ToolDescriptor tool(String name, Map inSchema = [type: 'object'], Map outSchema = null) {
        return new ToolDescriptor(name: name, description: "desc of ${name}", inputSchema: inSchema, outputSchema: outSchema)
    }

    def 'toolsFingerprint is null for no tools so the cache key stays byte-identical'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']]]
        def a = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))

        expect:
        AgentDef.toolsFingerprint(null, null) == null
        AgentDef.toolsFingerprint([] as List<ToolDescriptor>, [:]) == null

        and: 'a tool-free agent keeps an IDENTICAL key (no spurious resume invalidation)'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null) ==
            a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi')
    }

    def 'toolsFingerprint tracks the backing process script, not just the tool name'() {
        given:
        def upper = tool('uppercase')
        def v1 = AgentDef.toolsFingerprint([upper], [uppercase: 'tr a-z A-Z'])

        expect: 'editing the tool process script invalidates the agent cache entry'
        AgentDef.toolsFingerprint([upper], [uppercase: 'tr A-Z a-z']) != v1

        and: 'a changed descriptor (what the LLM is told about the tool) also invalidates it'
        AgentDef.toolsFingerprint([tool('uppercase', [type: 'object', properties: [text: [type: 'string']]])],
            [uppercase: 'tr a-z A-Z']) != v1

        and: 'an unchanged tool keeps the same fingerprint'
        AgentDef.toolsFingerprint([tool('uppercase')], [uppercase: 'tr a-z A-Z']) == v1
    }

    def 'toolsFingerprint is order-independent'() {
        given:
        def a = tool('alpha')
        def b = tool('beta')
        def sources = [alpha: 'body-a', beta: 'body-b']

        expect: 'the same set of tools, declared in either order, is the same key'
        AgentDef.toolsFingerprint([a, b], sources) == AgentDef.toolsFingerprint([b, a], sources)

        and: 'but a different set is a different key'
        AgentDef.toolsFingerprint([a], sources) != AgentDef.toolsFingerprint([a, b], sources)
    }

    def 'a sourceless tool (e.g. filesystem) still fingerprints its descriptor'() {
        expect: 'no backing process means a null source element, not a crash'
        AgentDef.toolsFingerprint([tool('filesystem')], [:]) != null
        AgentDef.toolsFingerprint([tool('filesystem')], [:]) == AgentDef.toolsFingerprint([tool('filesystem')], null)
    }

    // -----------------------------------------------------------------------
    // M-Endpoint resume: the RESOLVED endpoint folded into the cache key (design D5). A different
    // endpoint serves a different model under the same id, so a replay must not be shared across
    // endpoints -- while the credential is not part of an agent's identity and never enters.
    // -----------------------------------------------------------------------

    def 'canonicalAgentSource 7-arg with no endpoint equals the 6-arg form byte-for-byte'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']], required: ['x'], additionalProperties: false]
        def a = agentWith([instruction: 'be careful', goal: 'do it'], new PromptDef({ -> 'p' }, 'the-prompt'))
        def six = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null)

        expect: 'an agent with no endpoint set keeps an IDENTICAL key -- no spurious resume invalidation'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null) == six
        and: 'an empty endpoint is treated as unset (the core resolvers normalize "" to null anyway)'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, '') == six
        and: 'the byte-identical chain reaches the 3-arg form that predates all of this'
        six == a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi')
        !six.contains('baseUrl=')
    }

    def 'canonicalAgentSource folds in the resolved endpoint'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']]]
        def a = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))
        def none = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null)
        def local = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, 'http://localhost:8000/v1')

        expect: 'declaring an endpoint changes the key, appended as a trailing line'
        local != none
        local == none + '\nbaseUrl=http://localhost:8000/v1'

        and: 'a DIFFERENT endpoint is a different key'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, 'http://localhost:9000/v1') != local

        and: 'the endpoint composes with the tools fingerprint rather than replacing it'
        def withTools = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', 'tool-fp', 'http://localhost:8000/v1')
        withTools.contains('\ntools=tool-fp')
        withTools.endsWith('\nbaseUrl=http://localhost:8000/v1')
    }

    def 'no canonicalAgentSource overload can fold in a credential'() {
        given: 'drift guard for design D5 -- hashing a credential would invalidate every stored'
        // entry on key rotation, and would write a secret into the task hash inputs
        def overloads = AgentDef.declaredMethods.findAll { it.name == 'canonicalAgentSource' }

        expect: 'the widest overload is the 8-arg provider-aware one'
        overloads.size() == 6
        overloads*.parameterCount.max() == 8

        and: 'and its two trailing parameters are the endpoint and the credential NAMESPACE --'
        // the only two values AgentDef passes from the resolved settings, neither of them a secret
        overloads.find { it.parameterCount == 8 }.parameterTypes.toList() ==
            [String, Integer.TYPE, Map, List, String, String, String, String]
    }

    // -----------------------------------------------------------------------
    // D6: an EXPLICIT `agent.apiProvider` selects which environment variables the endpoint and the
    // credential come from, so it is part of how this agent was configured. An INFERRED one is a
    // pure function of `baseUrl`, which is already in the key -- and leaving it out means a later
    // addition to the D3 host table cannot silently invalidate anyone's stored runs.
    // -----------------------------------------------------------------------

    def 'canonicalAgentSource 8-arg with no apiProvider equals the 7-arg form byte-for-byte'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']], required: ['x'], additionalProperties: false]
        def a = agentWith([instruction: 'be careful', goal: 'do it'], new PromptDef({ -> 'p' }, 'the-prompt'))
        def seven = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null)

        expect: 'an agent that does not set the option keeps an IDENTICAL key -- no spurious invalidation'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, null) == seven
        and: 'an empty value is treated as unset (AgentConfig normalizes "" to null anyway)'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, '') == seven
        and: 'the byte-identical chain still reaches the narrowest form that carries a runner'
        seven == a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi')
        !seven.contains('apiProvider=')
        !seven.contains('baseUrl=')
    }

    def 'canonicalAgentSource folds in an explicit apiProvider, after the endpoint'() {
        given:
        def schema = [type: 'object', properties: [x: [type: 'string']]]
        def a = agentWith([instruction: 'i'], new PromptDef({ -> 'p' }, 'src'))
        def none = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, null)
        def routed = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, 'openrouter')

        expect: 'the option is appended as a trailing line, so it changes the key exactly once'
        routed == none + '\napiProvider=openrouter'

        and: 'a DIFFERENT namespace is a different key -- it selects a different credential'
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, 'azure') != routed

        and: 'setting it REDUNDANTLY to what the prefix already implied still changes the key'
        // the documented one-time invalidation: the key records the CONFIG, not the resolution
        a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', null, null, 'openai') != none

        and: 'it composes with the endpoint rather than replacing it, and comes after it'
        def both = a.canonicalAgentSource('openai/gpt-4o', 20, schema, null, 'pi', 'tool-fp', 'https://gw.corp/v1', 'openai')
        both.contains('\ntools=tool-fp')
        both.contains('\nbaseUrl=https://gw.corp/v1')
        both.endsWith('\napiProvider=openai')
    }

    def 'canonicalJson is key-sorted and insertion-order-independent'() {
        given:
        def m1 = [b: 1, a: 2, nested: [z: 1, y: 2]]
        def m2 = [a: 2, b: 1, nested: [y: 2, z: 1]]

        expect: 'reordered maps (including nested) produce the identical fingerprint'
        AgentDef.canonicalJson(m1) == AgentDef.canonicalJson(m2)

        and: 'keys are sorted in the output (a before b, y before z)'
        AgentDef.canonicalJson(m1) == '{"a":2,"b":1,"nested":{"y":2,"z":1}}'
    }

    def 'PromptDef 3-arg ctor stores valRefs; 2-arg ctor delegates to empty'() {
        given:
        def refs = [new TokenValRef('params.threshold'), new TokenValRef('task.ext.args')]

        expect: 'the 3-arg ctor carries valRefs and exposes their names'
        def p3 = new PromptDef({ -> 'x' }, 'src', refs)
        p3.valRefs == refs
        (p3.getValNames() as Set) == ['params.threshold', 'task.ext.args'] as Set

        and: 'the 2-arg ctor still works with empty valRefs'
        def p2 = new PromptDef({ -> 'x' }, 'src')
        p2.valRefs == []
        p2.getValNames() == []
    }

    def 'should expose the skills directive (single, list, none)'() {
        expect:
        new AgentDef(Mock(BaseScript), 'a', [skills: 'greet'] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == ['greet']
        new AgentDef(Mock(BaseScript), 'a', [skills: ['a', 'b']] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == ['a', 'b']
        new AgentDef(Mock(BaseScript), 'a', [:] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == []
    }

    def 'should no longer reject skills combined with a record (structured) output (M5 guard removed)'() {
        given:
        // M5 deletes the tools/skills-XOR-structured guard: skills + a record output is
        // now allowed (the plugin runs a final structuring turn). Proof at this unit level
        // is that run() proceeds PAST the (removed) guard to skills resolution instead of
        // short-circuiting with the old guard message.
        def inp = new AgentBuilder.AgentInput('q', String)
        def out = new AgentBuilder.AgentOutput('a', TestRec)
        def agent = new AgentDef(Mock(BaseScript), 'a', [skills: 'greet'] as Map<String,Object>, [inp], [out], new PromptDef({ -> 'h' }, 'h'))
        // inject a stub runner so run() gets PAST the runner lookup and deterministically
        // reaches skills resolution (otherwise it would abort earlier for a missing runner)
        nextflow.agent.AgentRunnerProvider.testRunner = { req -> '{}' } as nextflow.agent.AgentRunner

        when:
        agent.run(['x'] as Object[])

        then:
        // the old guard would have thrown ScriptRuntimeException with 'combining tools or
        // skills' synchronously; it is gone, so run() proceeds PAST it into skills resolution
        // and fails there because the `greet` skill is not on disk in this unit context. The
        // positive skills+structured record BIND is covered end-to-end in
        // AgentToolBridgeIntegrationTest 'should allow skills with a record (structured)
        // output and bind the JSON'.
        def e = thrown(ScriptRuntimeException)
        !(e.message?.contains('combining tools or skills'))
        and:
        // proves we reached skills resolution (the removed guard sat before it)
        e.message?.contains('skill') && e.message?.contains('greet')

        cleanup:
        nextflow.agent.AgentRunnerProvider.testRunner = null
    }
}
