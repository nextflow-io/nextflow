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

package nextflow.module

import java.nio.file.Path

import nextflow.exception.AbortOperationException
import org.yaml.snakeyaml.Yaml
import spock.lang.Specification
import spock.lang.TempDir

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleSpecFactoryTest extends Specification {

    @TempDir
    Path tempDir

    def 'should load valid spec' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''\
            name: nf-core/fastqc
            version: 1.0.0
            description: FastQC quality control
            keywords:
              - quality-control
              - fastq
            license: MIT
            authors:
              - John Doe
            requires:
              nextflow: ">=24.04.0"
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromYaml(metaYaml)

        then:
        spec.name == 'nf-core/fastqc'
        spec.version == '1.0.0'
        spec.description == 'FastQC quality control'
        spec.authors == ['John Doe']
        spec.license == 'MIT'
        spec.keywords == ['quality-control', 'fastq']
        spec.requires == ['nextflow': '>=24.04.0']
    }

    def 'should fail to load non-existent spec' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')

        when:
        ModuleSpecFactory.fromYaml(metaYaml)

        then:
        thrown(AbortOperationException)
    }

    def 'should load spec with standard inputs and outputs'() {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''\
            name: nf-core/align
            version: 1.0.0
            description: Align reads
            license: MIT
            input:
              - - name: meta
                  type: map
                  description: Sample metadata
                - name: reads
                  type: file
                  description: Input reads
              - name: index
                type: directory
                description: Index directory
            output:
              - - name: meta
                  type: map
                  description: Sample metadata
                - name: bam
                  type: file
                  description: Output BAM file
            topics:
              - - type: string
                  description: Process name
                - type: string
                  description: Tool name
                - type: string
                  description: Tool version
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromYaml(metaYaml)

        then:
        spec.inputs.size() == 2
        spec.inputs[0].isTuple()
        spec.inputs[0].components.size() == 2
        spec.inputs[0].components[0].name == 'meta'
        spec.inputs[0].components[0].type == 'map'
        spec.inputs[0].components[0].description == 'Sample metadata'
        spec.inputs[0].components[1].name == 'reads'
        spec.inputs[0].components[1].type == 'file'
        spec.inputs[0].components[1].description == 'Input reads'
        spec.inputs[1].name == 'index'
        spec.inputs[1].type == 'directory'
        spec.inputs[1].description == 'Index directory'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].components.size() == 2
        spec.outputs[0].components[0].name == 'meta'
        spec.outputs[0].components[0].type == 'map'
        spec.outputs[0].components[0].description == 'Sample metadata'
        spec.outputs[0].components[1].name == 'bam'
        spec.outputs[0].components[1].type == 'file'
        spec.outputs[0].components[1].description == 'Output BAM file'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[0].description == 'Process name'
        spec.topics[0].components[1].type == 'string'
        spec.topics[0].components[1].description == 'Tool name'
        spec.topics[0].components[2].type == 'string'
        spec.topics[0].components[2].description == 'Tool version'
    }

    def 'should load spec with nf-core style inputs and outputs'() {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''\
            name: nf-core/align
            version: 1.0.0
            description: Align reads
            license: MIT
            input:
              - - meta:
                    type: map
                    description: Sample metadata
                - reads:
                    type: file
                    description: Input reads
              - index:
                  type: directory
                  description: Index directory
            output:
              bam:
              - - meta:
                    type: map
                    description: Sample metadata
                - "*.bam":
                    type: file
                    description: Output BAM file
            topics:
              versions:
              - - process:
                    type: string
                    description: Process name
                - tool:
                    type: string
                    description: Tool name
                - version:
                    type: string
                    description: Tool version
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromYaml(metaYaml)

        then:
        spec.inputs.size() == 2
        spec.inputs[0].isTuple()
        spec.inputs[0].components.size() == 2
        spec.inputs[0].components[0].name == 'meta'
        spec.inputs[0].components[0].type == 'map'
        spec.inputs[0].components[0].description == 'Sample metadata'
        spec.inputs[0].components[1].name == 'reads'
        spec.inputs[0].components[1].type == 'file'
        spec.inputs[0].components[1].description == 'Input reads'
        spec.inputs[1].name == 'index'
        spec.inputs[1].type == 'directory'
        spec.inputs[1].description == 'Index directory'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].components.size() == 2
        spec.outputs[0].components[0].name == 'meta'
        spec.outputs[0].components[0].type == 'map'
        spec.outputs[0].components[0].description == 'Sample metadata'
        spec.outputs[0].components[1].name == '*.bam'
        spec.outputs[0].components[1].type == 'file'
        spec.outputs[0].components[1].description == 'Output BAM file'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[0].description == 'Process name'
        spec.topics[0].components[1].type == 'string'
        spec.topics[0].components[1].description == 'Tool name'
        spec.topics[0].components[2].type == 'string'
        spec.topics[0].components[2].description == 'Tool version'
    }

    def 'should load spec with no inputs or outputs'() {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''\
            name: nf-core/fastqc
            version: 1.0.0
            description: FastQC quality control
            license: MIT
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromYaml(metaYaml)

        then:
        !spec.inputs
        !spec.outputs
    }

    def 'should extract module spec from process definition'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process FASTQC {
                input:
                val sample_id
                path fasta
                tuple val(meta), path(reads)

                output:
                path "*.html"
                path "salmon" , emit: index
                tuple val(meta), path("*.bam")
                tuple val("${task.process}"), path('$prefix'), eval('salmon --version | sed -e "s/salmon //g"'), topic: versions

                script:
                "fastqc $reads"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.name == 'my-namespace/fastqc'

        and:
        spec.inputs.size() == 3
        spec.inputs[0].name == 'sample_id'
        spec.inputs[0].type == null
        spec.inputs[1].name == 'fasta'
        spec.inputs[1].type == 'file'
        spec.inputs[2].isTuple()
        spec.inputs[2].components.size() == 2
        spec.inputs[2].components[0].name == 'meta'
        spec.inputs[2].components[0].type == 'map'
        spec.inputs[2].components[1].name == 'reads'
        spec.inputs[2].components[1].type == 'file'

        and:
        spec.outputs.size() == 3
        spec.outputs[0].name == null
        spec.outputs[0].type == 'file'
        spec.outputs[1].name == 'index'
        spec.outputs[1].type == 'file'
        spec.outputs[2].isTuple()
        spec.outputs[2].components.size() == 2
        spec.outputs[2].components[0].name == 'meta'
        spec.outputs[2].components[0].type == 'map'
        spec.outputs[2].components[1].name == null
        spec.outputs[2].components[1].type == 'file'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].name == null
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[1].name == null
        spec.topics[0].components[1].type == 'file'
        spec.topics[0].components[2].name == null
        spec.topics[0].components[2].type == 'string'
    }

    def 'should omit input section when process has no input block'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process NO_INPUT {
                output:
                path "result.txt"

                script:
                "echo hello > result.txt"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.isEmpty()
        spec.outputs.size() == 1
    }

    def 'should extract and render module spec with missing fields'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process FASTQC {
                input:
                tuple val(meta), path(reads)

                output:
                path "*.html"

                script:
                "fastqc $reads"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')
        def yaml = spec.toYaml()
        def parsed = new Yaml().load(yaml) as Map

        then:
        // Top-level structure
        yaml.startsWith('# This file was auto-generated')
        parsed['name'] == 'my-namespace/fastqc'
        // null fields are omitted from the YAML map; TODO list is in the comment header
        !parsed.containsKey('version')
        yaml.contains('# TODO:')
        yaml.contains('Missing required field: version')
        yaml.contains('Missing required field: description')
        yaml.contains('Missing required field: license')

        and:
        // Tuple input is rendered as a nested list
        def input = parsed['input'] as List
        input.size() == 1
        input[0] instanceof List
        def tupleItems = input[0] as List
        tupleItems.size() == 2
        tupleItems[0] == [name: 'meta', type: 'map', description: 'TODO: Add description']
        tupleItems[1] == [name: 'reads', type: 'file', description: 'TODO: Add description']

        and:
        // Non-tuple output is a flat list entry
        def output = parsed['output'] as List
        output.size() == 1
        output[0] == [type: 'file', description: 'TODO: Add description']
    }

    def 'should throw error when no process is defined'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            workflow {
                channel.of(1,2,3).view()
            }
            '''.stripIndent()

        when:
        ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        thrown(AbortOperationException)
    }

    def 'should throw error when multiple processes are defined'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process FIRST {
                input:
                val x

                script:
                "echo $x"
            }

            process SECOND {
                input:
                path y

                script:
                "cat $y"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        thrown(AbortOperationException)
    }

    // =========================================================================
    // typed process tests
    // =========================================================================

    def 'should extract typed inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.enable.types = true

            process ALIGN {
                input:
                sample_id: String
                reads: List<Path>
                index: Path
                iteration: Integer

                output:
                bam: Path

                exec:
                bam = file("${sample_id}.bam")
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 4
        spec.inputs[0].name == 'sample_id'
        spec.inputs[0].type == 'string'
        spec.inputs[1].name == 'reads'
        spec.inputs[1].type == 'file'
        spec.inputs[2].name == 'index'
        spec.inputs[2].type == 'file'
        spec.inputs[3].name == 'iteration'
        spec.inputs[3].type == 'integer'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].name == 'bam'
        spec.outputs[0].type == 'file'
    }

    def 'should extract record inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.enable.types = true

            process ALIGN {
                input:
                record(sample_id: String, reads: List<Path>)

                output:
                record(sample_id: sample_id, bam: files("*.bam"))

                topic:
                record(process: "${task.process}", tool: 'align', version: '1.0.0') >> 'versions'

                exec:
                println "id=${sample_id}"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].isTuple()
        spec.inputs[0].components.size() == 2
        spec.inputs[0].components[0].name == 'sample_id'
        spec.inputs[0].components[0].type == 'string'
        spec.inputs[0].components[1].name == 'reads'
        spec.inputs[0].components[1].type == 'file'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].components.size() == 2
        spec.outputs[0].components[0].name == 'sample_id'
        spec.outputs[0].components[0].type == 'string'
        spec.outputs[0].components[1].name == 'bam'
        spec.outputs[0].components[1].type == 'file'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].name == 'process'
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[1].name == 'tool'
        spec.topics[0].components[1].type == 'string'
        spec.topics[0].components[2].name == 'version'
        spec.topics[0].components[2].type == 'string'
    }

    def 'should extract tuple inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.enable.types = true

            process ALIGN {
                input:
                tuple(sample_id: String, reads: List<Path>)

                output:
                tuple(sample_id, files("*.bam"))

                topic:
                tuple("${task.process}", 'align', '1.0.0') >> 'versions'

                exec:
                println "id=${sample_id}"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].isTuple()
        spec.inputs[0].components.size() == 2
        spec.inputs[0].components[0].name == 'sample_id'
        spec.inputs[0].components[0].type == 'string'
        spec.inputs[0].components[1].name == 'reads'
        spec.inputs[0].components[1].type == 'file'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].components.size() == 2
        spec.outputs[0].components[0].name == 'sample_id'
        spec.outputs[0].components[0].type == 'string'
        spec.outputs[0].components[1].name == null
        spec.outputs[0].components[1].type == 'file'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].name == null
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[1].name == null
        spec.topics[0].components[1].type == 'string'
        spec.topics[0].components[2].name == null
        spec.topics[0].components[2].type == 'string'
    }

    // =========================================================================
    // legacy process tests
    // =========================================================================

    def 'should extract legacy path inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process SAMTOOLS_SORT {
                input:
                path bam

                output:
                path "sorted.bam"

                script:
                "samtools sort $bam -o sorted.bam"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].name == 'bam'
        spec.inputs[0].type == 'file'
        and:
        spec.outputs.size() == 1
        spec.outputs[0].name == null
        spec.outputs[0].type == 'file'
    }

    def 'should extract legacy env inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process BWA_MEM {
                input:
                env 'GENOME_DIR'

                output:
                env 'ALIGNED_READS'

                script:
                "bwa mem \\$GENOME_DIR/ref.fa"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].name == null
        spec.inputs[0].type == 'string'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].name == null
        spec.outputs[0].type == 'string'
    }

    def 'should extract legacy stdin and stdout'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process CAT {
                input:
                stdin

                output:
                stdout

                script:
                "cat"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].name == null
        spec.inputs[0].type == 'string'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].name == null
        spec.outputs[0].type == 'string'
    }

    def 'should extract legacy val inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process COUNT {
                input:
                val message

                output:
                val count

                script:
                count = 42
                "echo hello"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].name == 'message'
        spec.inputs[0].type == null

        and:
        spec.outputs.size() == 1
        spec.outputs[0].name == 'count'
        spec.outputs[0].type == null
    }

    def 'should extract legacy tuple inputs and outputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process MIXED_TUPLE {
                input:
                tuple val(sample_id), path(reads), env('GENOME')

                output:
                tuple val(sample_id), path("*.bam"), env('ALIGNED_READS')
                tuple val("${task.process}"), val('align'), val('1.0.0'), topic: versions

                script:
                "bwa mem $reads"
            }
            '''.stripIndent()

        when:
        def spec = ModuleSpecFactory.fromScript(mainNf, namespace: 'my-namespace')

        then:
        spec.inputs.size() == 1
        spec.inputs[0].isTuple()
        spec.inputs[0].components.size() == 3
        spec.inputs[0].components[0].name == 'sample_id'
        spec.inputs[0].components[0].type == null
        spec.inputs[0].components[1].name == 'reads'
        spec.inputs[0].components[1].type == 'file'
        spec.inputs[0].components[2].name == null
        spec.inputs[0].components[2].type == 'string'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].components.size() == 3
        spec.outputs[0].components[0].name == 'sample_id'
        spec.outputs[0].components[0].type == null
        spec.outputs[0].components[1].name == null
        spec.outputs[0].components[1].type == 'file'
        spec.outputs[0].components[2].name == null
        spec.outputs[0].components[2].type == 'string'

        and:
        spec.topics.size() == 1
        spec.topics[0].isTuple()
        spec.topics[0].components.size() == 3
        spec.topics[0].components[0].name == null
        spec.topics[0].components[0].type == 'string'
        spec.topics[0].components[1].name == null
        spec.topics[0].components[1].type == 'string'
        spec.topics[0].components[2].name == null
        spec.topics[0].components[2].type == 'string'
    }

}
