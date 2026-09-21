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

import groovy.json.JsonOutput
import io.seqera.npr.api.schema.v1.ModuleChannel
import io.seqera.npr.api.schema.v1.ModuleChannelItem
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.api.schema.v1.ModuleTool
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam
import nextflow.script.AgentDef
import spock.lang.Specification

/**
 * <b>This is a tripwire, not a unit test.</b>
 *
 * <p>{@code AgentCacheKeyPinTest} pins {@code canonicalAgentSource}/{@code toolsFingerprint}, but
 * both are PURE FUNCTIONS OF THEIR ARGUMENTS, so a pin built from hand-made {@link ToolDescriptor}s
 * is green by construction no matter what happens to the code that PRODUCES those descriptors.
 * This file closes that hole. The chain it guards is:
 *
 * <pre>
 * ModuleSpecToolSchema.inputSchema / outputDescription  -&gt; ModuleToolBridge.wireSpec
 * ModuleMetadataToolSchema.inputSchema / description    -&gt; ToolDescriptor(name, description, inputSchema)
 *                                                      -&gt; AgentDef.toolsFingerprint hashes
 *                                                         d.description + canonicalJson(d.inputSchema)
 *                                                      -&gt; canonicalAgentSource `tools=` line
 *                                                      -&gt; BodyDef.source -&gt; task hash -&gt; -resume
 * </pre>
 *
 * <p>So the exact BYTES asserted below ARE the {@code -resume} cache key of every agent that wires a
 * module tool. A change to any of them silently invalidates every affected user's stored runs: the
 * pipeline still works, every other test still passes, and the only symptom is that resume re-runs
 * everything.
 *
 * <p><b>One hop in that chain is re-spelled, not executed.</b> {@code wireSpec} is an instance method
 * that needs a {@link nextflow.script.ProcessDef} (it reads the declared input-channel count off it),
 * which this file deliberately does not build -- the fixtures below are hand-made values, no session,
 * no process, no registry client. So the two schema producers and {@code buildDescription} ARE the
 * real code, reached directly or reflectively, but the {@code ToolDescriptor} they are assembled into
 * is composed HERE, by {@link #specDescriptor} / {@link #metadataDescriptor}, mirroring
 * {@code ModuleToolBridge.wireSpec}. A change to how {@code wireSpec} ITSELF composes the descriptor
 * -- a different tool name, a non-null {@code outputSchema}, a swapped source branch -- moves the
 * cache key and this file stays green. Keep the two in step by hand; if the descriptor ever becomes
 * buildable without a process, execute it here instead.
 *
 * <p>The expected values were obtained by RUNNING the code, not by reasoning about what it ought to
 * produce. They are asserted with {@code ==} against whole literals rather than {@code contains},
 * because the point is byte-identity — including the leading {@code 'Returns a JSON object with the
 * following output(s):'}, the trailing {@code 'File/path outputs are returned as absolute path
 * strings (never file contents).'}, the {@code '\n- `'} separators and the trailing space that an
 * empty tuple leaves behind.
 *
 * <p><b>Do not "fix" a failure here by updating an expected value.</b> A failure means a refactor
 * moved the agent task hash. Either revert the change, or — if the change is intentional and
 * accepted — treat it as a documented cache invalidation with a changelog entry, and say so.
 *
 * <p>Both input schemas are pinned twice: as {@code canonicalJson} (key-sorted — that is what
 * {@code toolsFingerprint} actually hashes) AND as plain {@code JsonOutput.toJson} (insertion order
 * — that is what the model is shown as the tool's {@code parameters}). The two envelopes and the
 * two type ladders deliberately DISAGREE between the {@code meta.yml} source and the registry
 * source (a {@code float} is {@code 'a string'} for one and {@code 'a number'} for the other; the
 * nf-core {@code meta.id} convention is unconditional for one and gated for the other), so both are
 * pinned in both variants.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleToolDescriptorPinTest extends Specification {

    // --- fixtures: hand-built values only, no session, no process, no registry client ---------

    static private ModuleParam p(String name, String type, String desc) {
        return new ModuleParam(name: name, type: type, description: desc)
    }

    static private ModuleParam tup(String name, ModuleParam... comps) {
        return new ModuleParam(name: name, components: comps.toList())
    }

    /**
     * A meta.yml-sourced module exercising every rung of the spec ladder: the unconditional
     * {@code meta} map convention, a file, a path, a directory, a plain map, an integer, a boolean,
     * a {@code val} and — the rung with no test of its own — a {@code float}, which this ladder
     * renders as a string. The outputs add a named tuple, an output with a NULL name, a tuple with
     * a NULL-named component and an EMPTY tuple (whose {@code 'an object with '} keeps a trailing
     * space).
     */
    static ModuleSpec richSpec() {
        return new ModuleSpec(
            name: 'RICH_MODULE',
            description: 'Assemble and quality-check a sample',
            inputs: [
                tup(null, p('meta', 'map', 'sample metadata'), p('reads', 'file', 'the reads')),
                p('reference', 'path', 'the reference genome'),
                p('threshold', 'float', 'the score threshold'),
                p('rounds', 'integer', null),
                p('strict', 'boolean', null),
                p('outdir', 'directory', 'where to write'),
                p('extras', 'map', 'free-form extras'),
                p('label', 'val', null),
            ],
            outputs: [
                tup('assembly', p('meta', 'map', null), p('fasta', 'path', 'the assembly')),
                p('score', 'float', 'the assembly score'),
                p('report', 'path', null),
                p(null, 'val', 'the unnamed one'),
                tup('mixed', p(null, 'file', null), p('n', 'integer', 'a count')),
                tup('empty'),
            ])
    }

    /** A module with NO inputs and NO outputs: pins the empty properties/required envelope. */
    static ModuleSpec bareSpec() {
        return new ModuleSpec(name: 'BARE_MODULE')
    }

    static private ModuleChannelItem item(String name, String type, String desc, String pattern = null, List<String> enumValues = null) {
        final result = new ModuleChannelItem().name(name).type(type).description(desc)
        if( pattern ) result.pattern(pattern)
        if( enumValues ) result._enum(enumValues)
        return result
    }

    static private ModuleChannel chan(boolean tuple, ModuleChannelItem... items) {
        return new ModuleChannel().tuple(tuple).items(items.toList())
    }

    /**
     * A registry-sourced module exercising every rung of the metadata ladder: the gated nf-core
     * {@code meta.id} convention, a file with a pattern, a path, an enum, an integer, a boolean, a
     * plain map and a {@code float}, which THIS ladder renders as a number. A null channel and a
     * null item pin the skips. The outputs add the {@code 'a value'} rung twice (a channel with no
     * items list and one with an empty items list), a multi-item tuple, a {@code tuple:true} emit
     * carrying a SINGLE item, a float item, a null-named item and a NULL emit name.
     */
    static ModuleMetadata richMetadata() {
        final out = new LinkedHashMap<String, ModuleChannel>()
        out.put('versions', new ModuleChannel())
        out.put('assembly', chan(true, item('meta', 'map', null), item('fasta', 'path', 'the assembly')))
        out.put('single', chan(true, item('bam', 'file', 'the alignment')))
        out.put('score', chan(false, item('score', 'float', 'the assembly score')))
        out.put('anon', chan(false, item(null, null, null)))
        out.put('empty', new ModuleChannel().items([]))
        out.put(null, chan(false, item('unnamed', 'string', 'the unnamed emit')))
        return new ModuleMetadata()
            .description('Assemble and quality-check a sample')
            .tools([
                new ModuleTool()
                    .name('spades')
                    .version('3.15.5')
                    .homepage(URI.create('https://example.org/spades'))
                    .documentation(URI.create('https://docs.example.org/spades')),
                new ModuleTool().name('quast'),
                new ModuleTool(),
            ])
            .input([
                chan(true, item('meta', 'map', 'sample metadata'), item('reads', 'file', 'the reads', '*.{fastq,fq}.gz')),
                chan(false, item('reference', 'path', 'the reference genome')),
                chan(false, item('threshold', 'float', 'the score threshold')),
                chan(false, item('rounds', 'integer', null)),
                chan(false, item('mode', 'string', 'the mode', null, ['fast', 'slow'])),
                chan(false, item('strict', 'boolean', null)),
                chan(false, item('extras', 'map', 'free-form extras')),
                null,
                new ModuleChannel().items([null]),
            ])
            .output(out)
    }

    /** Registry metadata with NO inputs, NO outputs, NO tools and NO description. */
    static ModuleMetadata bareMetadata() {
        return new ModuleMetadata()
    }

    // --- the two producers, assembled the way ModuleToolBridge.wireSpec assembles them. The
    //     producers are the real code; the assembly is a mirror -- see the class javadoc. --------

    /**
     * {@code AgentDef.canonicalJson} is {@code protected static} and {@code ModuleToolBridge
     * .buildDescription} is {@code private static}; both are reached reflectively so the pin runs
     * the REAL code rather than a re-spelling of it. A {@code NoSuchMethodException} here is itself
     * a finding: it means the descriptor-building path moved.
     */
    static private String canonicalJson(Object obj) {
        final m = AgentDef.getDeclaredMethod('canonicalJson', Object)
        m.accessible = true
        return (String) m.invoke(null, obj)
    }

    static private String specDescription(ModuleSpec spec) {
        final m = ModuleToolBridge.getDeclaredMethod('buildDescription', ModuleSpec)
        m.accessible = true
        return (String) m.invoke(null, spec)
    }

    /** Mirrors the descriptor {@code wireSpec} builds when there is no registry metadata. */
    static ToolDescriptor specDescriptor(String name, ModuleSpec spec) {
        return new ToolDescriptor(name, specDescription(spec), ModuleSpecToolSchema.inputSchema(spec), null)
    }

    /** Mirrors the descriptor {@code wireSpec} builds when the registry metadata IS the source. */
    static ToolDescriptor metadataDescriptor(String name, ModuleMetadata metadata, boolean nfCore) {
        return new ToolDescriptor(
            name,
            ModuleMetadataToolSchema.description(metadata),
            ModuleMetadataToolSchema.inputSchema(metadata, nfCore),
            null )
    }

    // --- the pins ----------------------------------------------------------------------------

    /**
     * An EMPTY tuple output renders as {@code 'an object with '} — {@code parts.join(', ')} over no
     * components leaves the trailing space behind. It is spelled as its own constant so no editor
     * or formatter can silently strip it out of the literal below.
     */
    static private final String EMPTY_TUPLE_LINE = '- `empty`: an object with '

    static private final String SPEC_OUTPUT_DESCRIPTION = '''\
Returns a JSON object with the following output(s):
- `assembly`: an object with `meta` (an object), `fasta` (a file path string) (the assembly)
- `score`: `score` (a string) (the assembly score)
- `report`: `report` (a file path string)
- `result`: `value` (a string) (the unnamed one)
- `mixed`: an object with `value` (a file path string), `n` (an integer) (a count)
''' + EMPTY_TUPLE_LINE + '''
File/path outputs are returned as absolute path strings (never file contents).'''

    static private final String META_OUTPUT_DESCRIPTION = '''\
Returns a JSON object with the following output(s):
- `versions`: a value
- `assembly`: an object with `meta` (an object), `fasta` (a file path string) (the assembly)
- `single`: an object with `bam` (a file path string) (the alignment)
- `score`: `score` (a number) (the assembly score)
- `anon`: `value` (a string)
- `empty`: a value
- `result`: `unnamed` (a string) (the unnamed emit)
File/path outputs are returned as absolute path strings (never file contents).'''

    static private final String META_TOOL_PREAMBLE = '''\
Assemble and quality-check a sample
Tool `spades` v3.15.5 (homepage: https://example.org/spades, documentation: https://docs.example.org/spades)
Tool `quast`'''

    def 'pins the meta.yml-sourced tool description, output prose included'() {
        given:
        def descriptor = specDescriptor('RICH_MODULE', richSpec())

        expect: 'the output prose stands alone byte-for-byte'
        ModuleSpecToolSchema.outputDescription(richSpec()) == SPEC_OUTPUT_DESCRIPTION

        and: 'and the descriptor is the module description, a blank line, then that prose'
        descriptor.name == 'RICH_MODULE'
        descriptor.description == 'Assemble and quality-check a sample\n\n' + SPEC_OUTPUT_DESCRIPTION
    }

    def 'pins the meta.yml-sourced input schema'() {
        given:
        def descriptor = specDescriptor('RICH_MODULE', richSpec())

        expect: 'key-sorted -- this is what toolsFingerprint hashes'
        canonicalJson(descriptor.inputSchema) == '{"additionalProperties":false,"properties":{"extras":{"additionalProperties":true,"description":"free-form extras","type":"object"},"label":{"type":"string"},"meta":{"additionalProperties":true,"description":"sample metadata","properties":{"id":{"description":"sample identifier","type":"string"}},"type":"object"},"outdir":{"description":"where to write (file path)","type":"string"},"reads":{"description":"the reads (file path)","type":"string"},"reference":{"description":"the reference genome (file path)","type":"string"},"rounds":{"type":"integer"},"strict":{"type":"boolean"},"threshold":{"description":"the score threshold","type":"string"}},"required":["meta","reads","reference","threshold","rounds","strict","outdir","extras","label"],"type":"object"}'

        and: 'insertion order -- this is what the model is shown as the tool parameters'
        JsonOutput.toJson(descriptor.inputSchema) == '{"type":"object","properties":{"meta":{"type":"object","description":"sample metadata","properties":{"id":{"type":"string","description":"sample identifier"}},"additionalProperties":true},"reads":{"type":"string","description":"the reads (file path)"},"reference":{"type":"string","description":"the reference genome (file path)"},"threshold":{"type":"string","description":"the score threshold"},"rounds":{"type":"integer"},"strict":{"type":"boolean"},"outdir":{"type":"string","description":"where to write (file path)"},"extras":{"type":"object","additionalProperties":true,"description":"free-form extras"},"label":{"type":"string"}},"required":["meta","reads","reference","threshold","rounds","strict","outdir","extras","label"],"additionalProperties":false}'
    }

    def 'pins the registry-sourced tool description, tools and output prose included'() {
        given:
        def descriptor = metadataDescriptor('RICH_MODULE', richMetadata(), true)

        expect: 'the output prose stands alone byte-for-byte'
        ModuleMetadataToolSchema.outputDescription(richMetadata()) == META_OUTPUT_DESCRIPTION

        and: 'and the descriptor is the module description, the tool lines, then that prose'
        descriptor.name == 'RICH_MODULE'
        descriptor.description == META_TOOL_PREAMBLE + '\n\n' + META_OUTPUT_DESCRIPTION

        and: 'the description does not depend on the nf-core flag'
        metadataDescriptor('RICH_MODULE', richMetadata(), false).description == descriptor.description
    }

    def 'pins the registry-sourced input schema, nf-core and not'() {
        given:
        def nfCore = metadataDescriptor('RICH_MODULE', richMetadata(), true)
        def plain = metadataDescriptor('RICH_MODULE', richMetadata(), false)

        expect: 'nf-core: `meta` carries the id convention'
        canonicalJson(nfCore.inputSchema) == '{"additionalProperties":false,"properties":{"extras":{"additionalProperties":true,"description":"free-form extras","type":"object"},"meta":{"additionalProperties":true,"description":"sample metadata","properties":{"id":{"description":"sample identifier","type":"string"}},"type":"object"},"mode":{"description":"the mode","enum":["fast","slow"],"type":"string"},"reads":{"description":"the reads (file path) (pattern: *.{fastq,fq}.gz)","type":"string"},"reference":{"description":"the reference genome (file path)","type":"string"},"rounds":{"type":"integer"},"strict":{"type":"boolean"},"threshold":{"description":"the score threshold","type":"number"}},"required":["meta","reads","reference","threshold","rounds","mode","strict","extras"],"type":"object"}'
        JsonOutput.toJson(nfCore.inputSchema) == '{"type":"object","properties":{"meta":{"type":"object","description":"sample metadata","properties":{"id":{"type":"string","description":"sample identifier"}},"additionalProperties":true},"reads":{"type":"string","description":"the reads (file path) (pattern: *.{fastq,fq}.gz)"},"reference":{"type":"string","description":"the reference genome (file path)"},"threshold":{"type":"number","description":"the score threshold"},"rounds":{"type":"integer"},"mode":{"type":"string","description":"the mode","enum":["fast","slow"]},"strict":{"type":"boolean"},"extras":{"type":"object","description":"free-form extras","additionalProperties":true}},"required":["meta","reads","reference","threshold","rounds","mode","strict","extras"],"additionalProperties":false}'

        and: 'not nf-core: `meta` is a plain open object'
        canonicalJson(plain.inputSchema) == '{"additionalProperties":false,"properties":{"extras":{"additionalProperties":true,"description":"free-form extras","type":"object"},"meta":{"additionalProperties":true,"description":"sample metadata","type":"object"},"mode":{"description":"the mode","enum":["fast","slow"],"type":"string"},"reads":{"description":"the reads (file path) (pattern: *.{fastq,fq}.gz)","type":"string"},"reference":{"description":"the reference genome (file path)","type":"string"},"rounds":{"type":"integer"},"strict":{"type":"boolean"},"threshold":{"description":"the score threshold","type":"number"}},"required":["meta","reads","reference","threshold","rounds","mode","strict","extras"],"type":"object"}'
        JsonOutput.toJson(plain.inputSchema) == '{"type":"object","properties":{"meta":{"type":"object","description":"sample metadata","additionalProperties":true},"reads":{"type":"string","description":"the reads (file path) (pattern: *.{fastq,fq}.gz)"},"reference":{"type":"string","description":"the reference genome (file path)"},"threshold":{"type":"number","description":"the score threshold"},"rounds":{"type":"integer"},"mode":{"type":"string","description":"the mode","enum":["fast","slow"]},"strict":{"type":"boolean"},"extras":{"type":"object","description":"free-form extras","additionalProperties":true}},"required":["meta","reads","reference","threshold","rounds","mode","strict","extras"],"additionalProperties":false}'
    }

    def 'pins the zero-input, zero-output envelope on both paths'() {
        given:
        def fromSpec = specDescriptor('BARE_MODULE', bareSpec())
        def fromMetadata = metadataDescriptor('BARE_MODULE', bareMetadata(), true)

        expect: 'an empty properties/required is EMITTED, never dropped'
        canonicalJson(fromSpec.inputSchema) == '{"additionalProperties":false,"properties":{},"required":[],"type":"object"}'
        JsonOutput.toJson(fromSpec.inputSchema) == '{"type":"object","properties":{},"required":[],"additionalProperties":false}'
        canonicalJson(fromMetadata.inputSchema) == '{"additionalProperties":false,"properties":{},"required":[],"type":"object"}'
        JsonOutput.toJson(fromMetadata.inputSchema) == '{"type":"object","properties":{},"required":[],"additionalProperties":false}'

        and: 'the meta.yml path falls back to the module name, the registry path to a literal'
        fromSpec.description == 'BARE_MODULE\n\nReturns a JSON object (no declared outputs).'
        fromMetadata.description == 'module tool\n\nReturns a JSON object (no declared outputs).'
    }
}
