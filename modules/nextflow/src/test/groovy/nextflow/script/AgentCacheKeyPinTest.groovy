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

import nextflow.agent.SkillDescriptor
import nextflow.agent.SkillResource
import nextflow.agent.ToolDescriptor
import spock.lang.Specification

/**
 * <b>This is a tripwire, not a unit test.</b>
 *
 * <p>{@link AgentDef#canonicalAgentSource} is written verbatim into {@code BodyDef.source}, which
 * {@code TaskHasher} folds into an agent task's hash; {@link AgentDef#toolsFingerprint} is one of
 * its lines. So the exact BYTES these two pure functions emit ARE the {@code -resume} cache key of
 * every agent ever run. A change to either silently invalidates every user's stored runs: the
 * pipeline still works, every other test still passes, and the only symptom is that resume re-runs
 * everything.
 *
 * <p>The expected values below were obtained by RUNNING the code, not by reasoning about what it
 * ought to produce. They are a record of what the parent commit does, deliberately asserted with
 * {@code ==} against whole literals rather than {@code contains}, because the point is
 * byte-identity — including the line order, the {@code \n} separators and the absence of a
 * trailing newline.
 *
 * <p><b>Do not "fix" a failure here by updating an expected value.</b> A failure means a refactor
 * changed the agent task hash. Either revert the change, or — if the change is intentional and
 * accepted — treat it as a documented cache invalidation with a changelog entry, and say so.
 *
 * <p>The fixture exercises the parts of the key that a refactor can plausibly disturb: a brokered
 * module tool with a backing process source, a second brokered tool with an empty input envelope
 * and no output schema and no source, a runner-native tool family with a runner identity, two
 * skills (one with a bundled resource), a wrapped multi-output schema whose keys are NOT in
 * alphabetical insertion order (so the canonical JSON sorting is pinned too), a resolved
 * {@code baseUrl} and an explicit {@code apiProvider}. Tools, skills and native refs are all
 * declared out of order, so their sorting is pinned as well.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentCacheKeyPinTest extends Specification {

    // --- the fixture: hand-built values only, no session, no runner, no dataflow -------------

    static private final String MODEL = 'anthropic/claude-sonnet-4-5-20250929'
    static private final int MAX_ITER = 20
    static private final String RUNNER = 'pi'
    static private final String RUNNER_ID = 'pi@0.5.0-alpha.1'
    static private final String BASE_URL = 'https://api.example.com/v1'
    static private final String API_PROVIDER = 'anthropic'
    static private final String INSTRUCTION = 'Be precise and cite the tool output.'
    static private final String GOAL = 'Produce a QC report for the sample.'
    static private final String PROMPT_SOURCE = 'Analyse ${sample} and report'

    /** A wrapper (multi-output) schema, keys deliberately NOT in alphabetical insertion order. */
    static private Map wrappedSchema() {
        return [
            type: 'object',
            properties: [
                summary: [type: 'string', description: 'A short summary'],
                report: [
                    type: 'object',
                    properties: [
                        score: [type: 'integer'],
                        label: [type: 'string', description: 'The label'] ],
                    required: ['score', 'label'],
                    additionalProperties: false ] ],
            required: ['summary', 'report'],
            additionalProperties: false ]
    }

    /** Two brokered module tools, declared out of name order. */
    static private List<ToolDescriptor> brokeredTools() {
        return [
            new ToolDescriptor(
                name: 'fastqc',
                description: 'Run FastQC on the reads.\nReturns a JSON object with the following output(s):\n- `html`: a file path string\nFile/path outputs are returned as absolute path strings (never file contents).',
                inputSchema: [
                    type: 'object',
                    properties: [
                        reads: [type: 'string', description: 'The reads to analyse'],
                        meta: [
                            type: 'object',
                            description: 'Sample metadata',
                            properties: [id: [type: 'string']],
                            additionalProperties: false ] ],
                    required: ['meta', 'reads'],
                    additionalProperties: false ],
                outputSchema: [
                    type: 'object',
                    properties: [html: [type: 'string']],
                    required: ['html'],
                    additionalProperties: false ]),
            new ToolDescriptor(
                name: 'align',
                description: 'Align the reads to the reference genome.',
                inputSchema: [
                    type: 'object',
                    properties: [:],
                    required: [],
                    additionalProperties: false ],
                outputSchema: null) ]
    }

    /** The backing {@code BodyDef.source} of each brokered tool; `align` deliberately has none. */
    static private Map<String,String> toolSources() {
        return [
            fastqc: '\n    fastqc --outdir . ${reads}\n    ',
            align: null ] as Map<String,String>
    }

    /** A runner-native tool family, declared out of order. */
    static private List<String> nativeRefs() {
        return ['fs:read', 'fs:grep']
    }

    /** Two skills, declared out of name order; one carries a bundled resource. */
    static private List<SkillDescriptor> skills() {
        return [
            new SkillDescriptor(
                name: 'variant-calling',
                description: 'How to call variants on this reference',
                content: '# Variant calling\n\nUse the joint caller.\n',
                resources: [new SkillResource('references/guide.md', '# Guide\n\nRead this first.\n')]),
            new SkillDescriptor(
                name: 'assembly',
                description: 'How to assemble a genome',
                content: '# Assembly\n',
                resources: null) ]
    }

    private AgentDef agent() {
        return new AgentDef(
            Mock(BaseScript),
            'qc_agent',
            [instruction: INSTRUCTION, goal: GOAL] as Map<String,Object>,
            [],
            [],
            new PromptDef({ -> 'p' }, PROMPT_SOURCE))
    }

    /** The digests below are the ones the code emits; they are quoted in the pinned strings too. */
    static private final String TOOLS_DIGEST = '2d816ebf596efd499f96b0ec3ac6a69e'
    static private final String SKILLS_DIGEST = '440b171d053115341f30ba92d4993717'

    // --- the pins ----------------------------------------------------------------------------

    def 'toolsFingerprint pins the digest of the brokered module tools'() {
        expect: 'the two descriptors are hashed name-sorted, with their backing process source'
        AgentDef.toolsFingerprint(brokeredTools(), toolSources()) == '42b9658507fd5d9ae72825c3e5cb2672'

        and: 'an agent with no tools at all still contributes nothing'
        AgentDef.toolsFingerprint(null, null) == null
    }

    def 'toolsFingerprint pins the digest of the runner-native tools'() {
        expect: 'the refs are hashed sorted, together with the runner identity, never a schema'
        AgentDef.toolsFingerprint(null, null, nativeRefs(), RUNNER_ID) == '7b3ed3b3ebdea669b70c0fb81126104a'
    }

    def 'toolsFingerprint pins the digest of the brokered and runner-native tools together'() {
        expect:
        AgentDef.toolsFingerprint(brokeredTools(), toolSources(), nativeRefs(), RUNNER_ID) == TOOLS_DIGEST
    }

    def 'skillsFingerprint pins the digest of the declared skills'() {
        expect: 'the descriptors are hashed name-sorted, over name/description/content/resources'
        AgentDef.skillsFingerprint(skills()) == SKILLS_DIGEST
    }

    def 'canonicalAgentSource pins the identity string of a tool-free, skill-free agent'() {
        given:
        def agent = agent()

        expect: 'byte-for-byte, no trailing newline, schema JSON key-sorted at every depth'
        agent.canonicalAgentSource(MODEL, MAX_ITER, wrappedSchema()) == '''\
agentModel=anthropic/claude-sonnet-4-5-20250929
temperature=default
instruction=Be precise and cite the tool output.
goal=Produce a QC report for the sample.
maxIterations=20
prompt=Analyse ${sample} and report
outputSchema={"additionalProperties":false,"properties":{"report":{"additionalProperties":false,"properties":{"label":{"description":"The label","type":"string"},"score":{"type":"integer"}},"required":["score","label"],"type":"object"},"summary":{"description":"A short summary","type":"string"}},"required":["summary","report"],"type":"object"}'''
    }

    def 'canonicalAgentSource pins the identity string of an agent with a runner, tools, skills, an endpoint and a provider'() {
        given:
        def agent = agent()
        def tools = AgentDef.toolsFingerprint(brokeredTools(), toolSources(), nativeRefs(), RUNNER_ID)

        expect: 'byte-for-byte: the leading agentRunner line, then the four trailing lines in this order'
        agent.canonicalAgentSource(MODEL, MAX_ITER, wrappedSchema(), skills(), RUNNER, tools, BASE_URL, API_PROVIDER) == '''\
agentRunner=pi
agentModel=anthropic/claude-sonnet-4-5-20250929
temperature=default
instruction=Be precise and cite the tool output.
goal=Produce a QC report for the sample.
maxIterations=20
prompt=Analyse ${sample} and report
outputSchema={"additionalProperties":false,"properties":{"report":{"additionalProperties":false,"properties":{"label":{"description":"The label","type":"string"},"score":{"type":"integer"}},"required":["score","label"],"type":"object"},"summary":{"description":"A short summary","type":"string"}},"required":["summary","report"],"type":"object"}
skills=440b171d053115341f30ba92d4993717
tools=2d816ebf596efd499f96b0ec3ac6a69e
baseUrl=https://api.example.com/v1
apiProvider=anthropic'''
    }
}
