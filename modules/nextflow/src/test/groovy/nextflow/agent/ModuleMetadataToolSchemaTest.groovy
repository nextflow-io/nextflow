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

import io.seqera.npr.api.schema.v1.ModuleChannel
import io.seqera.npr.api.schema.v1.ModuleChannelItem
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.api.schema.v1.ModuleTool
import spock.lang.Specification

/**
 * Verifies that {@link ModuleMetadataToolSchema} maps the registry {@link ModuleMetadata}
 * into a portable JSON-schema (with per-field descriptions, patterns, enums and the nf-core
 * {@code meta.id} convention) and a human-readable tool description.
 *
 * The npr-api DTOs ARE test-constructible: each has a public no-arg constructor plus fluent
 * builders (e.g. {@code name(..)}/{@code type(..)}/{@code _enum(..)}/{@code items(..)}), so the
 * test builds a real {@link ModuleMetadata} rather than a stub.
 */
class ModuleMetadataToolSchemaTest extends Specification {

    private static ModuleChannelItem item(String name, String type, String desc, String pattern = null, List<String> enumValues = null) {
        def it = new ModuleChannelItem().name(name).type(type).description(desc)
        if( pattern ) it.pattern(pattern)
        if( enumValues ) it._enum(enumValues)
        return it
    }

    private static ModuleChannel tuple(ModuleChannelItem... items) {
        return new ModuleChannel().tuple(true).items(items.toList())
    }

    private static ModuleChannel scalar(ModuleChannelItem item) {
        return new ModuleChannel().tuple(false).items([item])
    }

    def 'should flatten an nf-core tuple input with the meta.id convention'() {
        given:
        def metadata = new ModuleMetadata()
            .description('Run FastQC on sequenced reads')
            .input([ tuple(
                item('meta', 'map', 'sample meta'),
                item('reads', 'file', 'input reads', '*.{fastq,fq}.gz') ) ])

        when:
        def schema = ModuleMetadataToolSchema.inputSchema(metadata, true)

        then:
        schema.type == 'object'
        schema.additionalProperties == false
        schema.required == ['meta', 'reads']
        and:
        // meta.id convention -> nested object with an id string property
        schema.properties.meta.type == 'object'
        schema.properties.meta.properties.id.type == 'string'
        schema.properties.meta.properties.id.description == 'sample identifier'
        schema.properties.meta.additionalProperties == true
        and:
        // file -> string path, with the pattern appended to the description
        schema.properties.reads.type == 'string'
        schema.properties.reads.description == 'input reads (file path) (pattern: *.{fastq,fq}.gz)'
    }

    def 'should NOT apply the meta.id convention for a non-nf-core module'() {
        given:
        def metadata = new ModuleMetadata()
            .input([ scalar(item('meta', 'map', 'sample meta')) ])

        when:
        def schema = ModuleMetadataToolSchema.inputSchema(metadata, false)

        then:
        // generic map -> open object, no nested id property
        schema.properties.meta.type == 'object'
        schema.properties.meta.additionalProperties == true
        schema.properties.meta.description == 'sample meta'
        !schema.properties.meta.containsKey('properties')
    }

    def 'should map scalar / integer / number / boolean / enum input types'() {
        given:
        def metadata = new ModuleMetadata().input([
            scalar(item('label', 'string', 'a label')),
            scalar(item('count', 'integer', 'a count')),
            scalar(item('ratio', 'float', 'a ratio')),
            scalar(item('flag', 'boolean', null)),
            scalar(item('mode', 'string', 'a mode', null, ['fast', 'slow'])),
        ])

        when:
        def schema = ModuleMetadataToolSchema.inputSchema(metadata, false)

        then:
        schema.properties.label.type == 'string'
        schema.properties.count.type == 'integer'
        schema.properties.ratio.type == 'number'
        schema.properties.flag.type == 'boolean'
        and:
        schema.properties.mode.type == 'string'
        schema.properties.mode.enum == ['fast', 'slow']
        and:
        schema.required == ['label', 'count', 'ratio', 'flag', 'mode']
    }

    def 'should throw on a duplicate flattened property name'() {
        given:
        def metadata = new ModuleMetadata().input([
            tuple(item('meta', 'map', null), item('reads', 'file', null)),
            scalar(item('meta', 'map', null)),
        ])

        when:
        ModuleMetadataToolSchema.inputSchema(metadata, true)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('duplicate input name `meta`')
    }

    def 'should throw on an input item with no name'() {
        given:
        def metadata = new ModuleMetadata().input([ scalar(item(null, 'file', 'no name')) ])

        when:
        ModuleMetadataToolSchema.inputSchema(metadata, false)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should build the description from module + tools + output shape'() {
        given:
        def metadata = new ModuleMetadata()
            .description('Run FastQC on sequenced reads')
            .tools([ new ModuleTool()
                        .name('fastqc')
                        .version('0.12.1')
                        .homepage(URI.create('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/')) ])
            .output([
                report: tuple(item('meta', 'map', null), item('html', 'file', 'the report html')),
            ])

        when:
        def text = ModuleMetadataToolSchema.description(metadata)

        then:
        text.contains('Run FastQC on sequenced reads')
        text.contains('fastqc')
        text.contains('0.12.1')
        text.contains('homepage: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/')
        and:
        text.contains('`report`')
        text.contains('`html`')
        text.contains('the report html')
        text.contains('a file path string')
        text.contains('absolute path strings')
    }

    def 'should report flattened input property names'() {
        given:
        def metadata = new ModuleMetadata().input([
            tuple(item('meta', 'map', null), item('reads', 'file', null)),
            scalar(item('db', 'file', null)),
        ])

        expect:
        ModuleMetadataToolSchema.inputPropertyNames(metadata) == ['meta', 'reads', 'db']
    }

    def 'should handle null metadata / empty inputs gracefully'() {
        expect:
        ModuleMetadataToolSchema.inputSchema(null, false).properties == [:]
        ModuleMetadataToolSchema.inputPropertyNames(null) == []
        ModuleMetadataToolSchema.description(null).contains('module tool')
        ModuleMetadataToolSchema.outputDescription(null).contains('no declared outputs')
    }
}
