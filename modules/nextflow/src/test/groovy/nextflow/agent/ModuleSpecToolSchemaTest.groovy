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

import java.nio.file.Files

import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam
import nextflow.module.ModuleSpecFactory
import spock.lang.Specification

/**
 * Verifies that {@link ModuleSpecToolSchema} flattens a module {@code meta.yml}
 * spec into a portable JSON-schema with per-component properties + descriptions,
 * and renders a human-readable output description.
 */
class ModuleSpecToolSchemaTest extends Specification {

    private static ModuleParam param(String name, String type, String desc) {
        return new ModuleParam(name: name, type: type, description: desc)
    }

    private static ModuleParam tuple(ModuleParam... comps) {
        return new ModuleParam(components: comps.toList())
    }

    def 'should flatten a tuple input channel into per-component properties with descriptions'() {
        given:
        def spec = new ModuleSpec(
            name: 'fastqc',
            description: 'Run FastQC',
            inputs: [ tuple(param('meta', 'map', 'sample meta'), param('reads', 'file', 'input reads')) ],
            outputs: [
                new ModuleParam(name: 'zip', components: [param('meta', 'map', null), param('zip', 'file', 'the zip')]),
            ]
        )

        when:
        def schema = ModuleSpecToolSchema.inputSchema(spec)

        then:
        schema.type == 'object'
        schema.additionalProperties == false
        schema.required == ['meta', 'reads']
        and:
        schema.properties.meta.type == 'object'
        schema.properties.meta.additionalProperties == true
        schema.properties.meta.description == 'sample meta'
        // nf-core meta.id convention: meta input carries a nested `id` property
        schema.properties.meta.properties.id.type == 'string'
        schema.properties.meta.properties.id.description == 'sample identifier'
        and:
        schema.properties.reads.type == 'string'
        schema.properties.reads.description == 'input reads (file path)'
    }

    def 'should apply nf-core meta.id convention for map input named meta'() {
        given:
        def spec = new ModuleSpec(
            name: 'fastqc',
            inputs: [ tuple(param('meta', 'map', 'sample meta'), param('reads', 'file', 'input reads')) ]
        )

        when:
        def schema = ModuleSpecToolSchema.inputSchema(spec)

        then:
        schema.properties.meta.type == 'object'
        schema.properties.meta.description == 'sample meta'
        schema.properties.meta.properties.id.type == 'string'
        schema.properties.meta.properties.id.description == 'sample identifier'
        schema.properties.meta.additionalProperties == true
        and:
        schema.properties.reads.type == 'string'
    }

    def 'should map scalar / integer / boolean input types leniently'() {
        given:
        def spec = new ModuleSpec(
            name: 'mod',
            inputs: [
                param('label', 'val', 'a label'),
                param('count', 'integer', 'a count'),
                param('flag', 'boolean', null),
                param('weird', 'something-unknown', 'fallback'),
            ]
        )

        when:
        def schema = ModuleSpecToolSchema.inputSchema(spec)

        then:
        schema.properties.label.type == 'string'
        schema.properties.count.type == 'integer'
        schema.properties.flag.type == 'boolean'
        schema.properties.weird.type == 'string'
        schema.required == ['label', 'count', 'flag', 'weird']
    }

    def 'should throw on duplicate flattened property names across channels'() {
        given:
        def spec = new ModuleSpec(
            name: 'mod',
            inputs: [
                tuple(param('meta', 'map', null), param('reads', 'file', null)),
                param('meta', 'map', null),
            ]
        )

        when:
        ModuleSpecToolSchema.inputSchema(spec)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('duplicate input name `meta`')
    }

    def 'should describe outputs as prose including file path contract'() {
        given:
        def spec = new ModuleSpec(
            name: 'fastqc',
            outputs: [
                new ModuleParam(name: 'report', components: [param('meta', 'map', null), param('outfile', 'file', 'the output')]),
            ]
        )

        when:
        def text = ModuleSpecToolSchema.outputDescription(spec)

        then:
        text.contains('`report`')
        text.contains('`outfile`')
        text.contains('the output')
        text.contains('a file path string')
        text.contains('absolute path strings')
    }

    def 'should derive the flattened schema from a fastqc-style meta.yml'() {
        given:
        def dir = Files.createTempDirectory('test')
        def meta = dir.resolve('meta.yml')
        meta.text = '''\
            name: fastqc
            description: Run FastQC on sequenced reads
            input:
              - - meta:
                    type: map
                    description: sample meta
                - reads:
                    type: file
                    description: input reads
            output:
              report:
                - - meta:
                      type: map
                      description: sample meta
                  - outfile:
                      type: file
                      description: the output report
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromYaml(meta)
        def schema = ModuleSpecToolSchema.inputSchema(spec)

        then:
        schema.properties.meta.type == 'object'
        schema.properties.meta.description == 'sample meta'
        schema.properties.reads.type == 'string'
        schema.properties.reads.description == 'input reads (file path)'
        schema.required == ['meta', 'reads']
        and:
        // NOTE: ModuleSpecFactory.fromYaml drops the Map output channel key (`report`),
        // so the channel-level emit name is not recovered from a Map-keyed meta.yml; the
        // tuple component names (meta/outfile) are preserved. The bridge recovers the real
        // emit name from the process ChannelOut at runtime.
        def desc = ModuleSpecToolSchema.outputDescription(spec)
        desc.contains('`meta`')
        desc.contains('`outfile`')
        desc.contains('the output report')
        desc.contains('absolute path strings')
    }
}
