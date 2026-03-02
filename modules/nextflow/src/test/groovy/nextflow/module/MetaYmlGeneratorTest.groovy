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

import nextflow.exception.AbortOperationException
import org.yaml.snakeyaml.Yaml
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Tests for MetaYmlGenerator
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class MetaYmlGeneratorTest extends Specification {

    @TempDir
    Path tempDir

    // =========================================================================
    // V1 classic qualifier syntax tests (T008)
    // =========================================================================

    def 'should extract V1 process with val and path inputs'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process FASTQC {
                input:
                val sample_id
                path reads

                output:
                path "*.html"

                script:
                "fastqc $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.name == 'FASTQC'
        meta.syntaxVersion == 'V1'
        meta.inputs.size() == 2
        meta.inputs[0].identifier == 'sample_id'
        meta.inputs[0].resolvedType == 'string'
        meta.inputs[1].identifier == 'reads'
        meta.inputs[1].resolvedType == 'file'
        meta.outputs.size() == 1
        meta.outputs[0].identifier == '*.html'
        meta.outputs[0].resolvedType == 'list'
    }

    def 'should extract V1 tuple input with val(meta) and path(reads)'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ALIGN {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("*.bam")

                script:
                "bwa $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.size() == 1
        meta.inputs[0].isTuple()
        meta.inputs[0].subEntries.size() == 2
        meta.inputs[0].subEntries[0].identifier == 'meta'
        meta.inputs[0].subEntries[0].resolvedType == 'map'
        meta.inputs[0].subEntries[1].identifier == 'reads'
        meta.inputs[0].subEntries[1].resolvedType == 'file'

        and:
        meta.outputs.size() == 1
        meta.outputs[0].isTuple()
        meta.outputs[0].subEntries.size() == 2
        meta.outputs[0].subEntries[0].identifier == 'meta'
        meta.outputs[0].subEntries[0].resolvedType == 'map'
        meta.outputs[0].subEntries[1].identifier == '*.bam'
        meta.outputs[0].subEntries[1].resolvedType == 'list'
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
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.isEmpty()
        meta.outputs.size() == 1
    }

    def 'should lowercase process name in generated name field'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process SAMTOOLS_SORT {
                input:
                path bam

                script:
                "samtools sort $bam"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)
        def yaml = MetaYmlGenerator.render(meta)

        then:
        meta.name == 'SAMTOOLS_SORT'
        yaml.contains('<scope-placeholder>/samtools_sort')
    }

    def 'should render V1 YAML with correct structure'() {
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
        def meta = MetaYmlGenerator.extract(mainNf)
        def yaml = MetaYmlGenerator.render(meta)
        def parsed = new Yaml().load(yaml) as Map

        then:
        // Top-level structure
        yaml.startsWith('# This file was auto-generated')
        parsed['name'] == '<scope-placeholder>/fastqc'
        parsed.containsKey('version')

        and:
        // Tuple input is rendered as a nested list
        def input = parsed['input'] as List
        input.size() == 1
        input[0] instanceof List
        def tupleItems = input[0] as List
        tupleItems.size() == 2
        tupleItems[0] == [meta: [type: 'map', description: 'TODO: Add description']]
        tupleItems[1] == [reads: [type: 'file', description: 'TODO: Add description']]

        and:
        // Non-tuple output is a flat list entry
        def output = parsed['output'] as List
        output.size() == 1
        output[0] == ['*.html': [type: 'list', description: 'TODO: Add description']]
    }

    def 'should warn and use first process when multiple are defined'() {
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
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.name == 'FIRST'
    }

    def 'should throw AbortOperationException when no process found'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            workflow {
                Channel.of(1,2,3) | view
            }
            '''.stripIndent()

        when:
        MetaYmlGenerator.extract(mainNf)

        then:
        thrown(AbortOperationException)
    }

    // =========================================================================
    // V2 static type syntax tests (T009)
    // =========================================================================

    def 'should extract V2 typed inputs with simple types'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process ALIGN {
                input:
                sample_id: String
                reads: Path
                index: int

                output:
                bam: Path

                exec:
                bam = file("${sample_id}.bam")
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.syntaxVersion == 'V2'
        meta.inputs.size() == 3
        meta.inputs[0].identifier == 'sample_id'
        meta.inputs[0].resolvedType == 'string'
        meta.inputs[1].identifier == 'reads'
        meta.inputs[1].resolvedType == 'file'
        meta.inputs[2].identifier == 'index'
        meta.inputs[2].resolvedType == 'integer'
    }

    def 'should extract V2 output with typed declaration'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process ALIGN {
                input:
                reads: Path

                output:
                bam: Path

                exec:
                bam = file("out.bam")
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bam'
        meta.outputs[0].resolvedType == 'file'
    }

    def 'should extract V2 tuple input with TupleParameter'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process ALIGN {
                input:
                (sample_id, reads): Tuple<String, Path>
                index: Path
                id: String
                anon_tuple: Tuple<String, Path>

                output:
                bam = file("out.bam")
                out_tuple = tuple(sample_id, file("out2.bam"))
                res = test
                res2 = id

                exec:
                    test = "hola"

            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        // 3 inputs: tuple(sample_id, reads), index, id
        meta.inputs.size() == 4
        meta.inputs[0].isTuple()
        meta.inputs[0].subEntries.size() == 2
        meta.inputs[0].subEntries[0].identifier == 'sample_id'
        meta.inputs[0].subEntries[0].resolvedType == 'string'
        meta.inputs[0].subEntries[1].identifier == 'reads'
        meta.inputs[0].subEntries[1].resolvedType == 'file'
        meta.inputs[1].identifier == 'index'
        meta.inputs[1].resolvedType == 'file'
        meta.inputs[2].identifier == 'id'
        meta.inputs[2].resolvedType == 'string'
        !meta.inputs[3].isTuple()
        meta.inputs[3].identifier == 'anon_tuple'
        meta.inputs[3].resolvedType == 'tuple'

        and:
        // 4 outputs inferred via symbol map
        meta.outputs.size() == 4
        // bam = file("out.bam") → RHS is file() call → file
        meta.outputs[0].identifier == 'bam'
        meta.outputs[0].resolvedType == 'file'
        // out_tuple = tuple(sample_id, file("out2.bam")) → tuple with 2 sub-entries
        meta.outputs[1].identifier == 'out_tuple'
        meta.outputs[1].isTuple()
        meta.outputs[1].subEntries.size() == 2
        meta.outputs[1].subEntries[0].identifier == 'sample_id'
        meta.outputs[1].subEntries[0].resolvedType == 'string'   // from input symbol map
        meta.outputs[1].subEntries[1].identifier == 'out2.bam'  // from file("out2.bam") arg
        meta.outputs[1].subEntries[1].resolvedType == 'file'
        // res = test → exec: test = "hola" → symbol map: test → string
        meta.outputs[2].identifier == 'res'
        meta.outputs[2].resolvedType == 'string'
        // res2 = id → input id: String → symbol map: id → string
        meta.outputs[3].identifier == 'res2'
        meta.outputs[3].resolvedType == 'string'
    }

    def 'should extract V2 output with assignment expression'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process MY_PROCESS {
                input:
                typed_input: Path
                id: String
                iteration: int

                output:
                result = tuple(id, file("bam.out"))
                outIteration = iteration

                exec:
                test = 'done'
                result = test
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs[0].identifier == 'typed_input'
        meta.inputs[0].resolvedType == 'file'
        // exec: result = test; test = 'done' → symbol map has result → string
        // so exec block assignment overrides the tuple() in the output block
        meta.outputs.size() == 2
        meta.outputs[0].identifier == 'result'
        meta.outputs[0].resolvedType == 'string'
        // outIteration = iteration → declared type String → string
        meta.outputs[1].identifier == 'outIteration'
        meta.outputs[1].resolvedType == 'integer'
    }

    // =========================================================================
    // V1 path() output type resolution tests (T010)
    // =========================================================================

    def 'should resolve V1 path output with literal filename as file'() {
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
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'sorted.bam'
        meta.outputs[0].resolvedType == 'file'
    }

    def 'should resolve V1 path output with glob pattern as list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = """\
            process P {
                output:
                path "${pattern}"

                script: "echo"
            }
            """.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].resolvedType == 'list'

        where:
        pattern << ['*.bam', 'result?.bam', '{reads,index}.fa', '[0-9]*.txt']
    }

    def 'should resolve V1 path output with arity 1 as file'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process SAMTOOLS_SORT {
                input:
                path bam

                output:
                path "*.bam", arity: '1'

                script:
                "samtools sort $bam -o sorted.bam"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].resolvedType == 'file'
    }

    def 'should resolve V1 path output with arity range as list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process SAMTOOLS_SORT {
                input:
                path bam

                output:
                path "output.bam", arity: '1..*'

                script:
                "samtools sort $bam"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].resolvedType == 'list'
    }

    def 'should resolve V1 path input variable as file regardless of name'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ALIGN {
                input:
                path reads

                script:
                "bwa $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.size() == 1
        meta.inputs[0].identifier == 'reads'
        meta.inputs[0].resolvedType == 'file'
    }

    // =========================================================================
    // V1 additional qualifier tests
    // =========================================================================

    def 'should extract V1 env input and env output'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ENV_PROCESS {
                input:
                env GENOME_DIR

                output:
                env ALIGNED_READS

                script:
                "bwa mem $GENOME_DIR/ref.fa"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.size() == 1
        meta.inputs[0].identifier == 'GENOME_DIR'
        meta.inputs[0].resolvedType == 'string'

        and:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'ALIGNED_READS'
        meta.outputs[0].resolvedType == 'string'
    }

    def 'should extract V1 stdin input'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process STDIN_PROCESS {
                input:
                stdin

                script:
                "cat"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.size() == 1
        meta.inputs[0].identifier == 'stdin'
        meta.inputs[0].resolvedType == 'stdin'
    }

    def 'should extract V1 stdout output'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ECHO_PROCESS {
                input:
                val message

                output:
                stdout

                script:
                "echo $message"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'stdout'
        meta.outputs[0].resolvedType == 'stdout'
    }

    def 'should extract V1 val output'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process COUNT_READS {
                input:
                path bam

                output:
                val count

                script:
                "echo hello"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'count'
        meta.outputs[0].resolvedType == 'string'
    }

    def 'should extract V1 deprecated file qualifier as file and list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process FILE_PROCESS {
                input:
                file reads

                output:
                file "result.txt"
                file "*.log"

                script:
                "process $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs.size() == 1
        meta.inputs[0].identifier == 'reads'
        meta.inputs[0].resolvedType == 'file'

        and:
        meta.outputs.size() == 2
        meta.outputs[0].identifier == 'result.txt'
        meta.outputs[0].resolvedType == 'file'
        meta.outputs[1].identifier == '*.log'
        meta.outputs[1].resolvedType == 'list'
    }

    def 'should extract V1 tuple with val, path and env sub-qualifiers'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process MIXED_TUPLE {
                input:
                tuple val(sample_id), path(reads), env(GENOME)

                output:
                tuple val(sample_id), path("*.bam"), env(ALIGNED_READS)

                script:
                "bwa mem $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.inputs[0].isTuple()
        meta.inputs[0].subEntries.size() == 3
        meta.inputs[0].subEntries[0].identifier == 'sample_id'
        meta.inputs[0].subEntries[0].resolvedType == 'string'
        meta.inputs[0].subEntries[1].identifier == 'reads'
        meta.inputs[0].subEntries[1].resolvedType == 'file'
        meta.inputs[0].subEntries[2].identifier == 'GENOME'
        meta.inputs[0].subEntries[2].resolvedType == 'string'

        and:
        meta.outputs[0].isTuple()
        meta.outputs[0].subEntries.size() == 3
        meta.outputs[0].subEntries[0].identifier == 'sample_id'
        meta.outputs[0].subEntries[0].resolvedType == 'string'
        meta.outputs[0].subEntries[1].identifier == '*.bam'
        meta.outputs[0].subEntries[1].resolvedType == 'list'
        meta.outputs[0].subEntries[2].identifier == 'ALIGNED_READS'
        meta.outputs[0].subEntries[2].resolvedType == 'string'
    }

    // =========================================================================
    // V2 path() output type resolution tests
    // =========================================================================

    def 'should resolve V2 path() output with glob as list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process SAMTOOLS {
                input:
                reads: Path

                output:
                bam = path("*.bam")

                script:
                "samtools sort $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bam'
        meta.outputs[0].resolvedType == 'list'
    }

    def 'should resolve V2 path() output with literal name as file'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process SAMTOOLS {
                input:
                reads: Path

                output:
                bam = path("sorted.bam")

                script:
                "samtools sort $reads -o sorted.bam"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bam'
        meta.outputs[0].resolvedType == 'file'
    }

    def 'should resolve V2 path() output with arity 1 as file'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process SAMTOOLS {
                input:
                reads: Path

                output:
                bam = path("*.bam", arity: '1')

                script:
                "samtools sort $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bam'
        meta.outputs[0].resolvedType == 'file'
    }

    def 'should resolve V2 path() output with arity range as list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process SAMTOOLS {
                input:
                reads: Path

                output:
                bams = path("*.bam", arity: '1..*')

                script:
                "samtools sort $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bams'
        meta.outputs[0].resolvedType == 'list'
    }

    def 'should resolve V2 env() output as string'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process ENV_PROCESS {
                input:
                reads: Path

                output:
                result = env(MY_VAR)

                script:
                "export MY_VAR=done"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'result'
        meta.outputs[0].resolvedType == 'string'
    }

    def 'should resolve V2 stdout() output as string'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            nextflow.preview.types = true

            process ECHO_PROCESS {
                input:
                msg: String

                output:
                result = stdout()

                script:
                "echo $msg"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'result'
        meta.outputs[0].resolvedType == 'string'
    }

    def 'should use emit name as tuple identifier in V1 output'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ALIGN {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("*.bam"), emit: bam_ch

                script:
                "bwa $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)

        then:
        meta.outputs.size() == 1
        meta.outputs[0].identifier == 'bam_ch'
        meta.outputs[0].isTuple()
        meta.outputs[0].subEntries[1].resolvedType == 'list'
    }

    def 'should render named V1 tuple output under its emit identifier'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ALIGN {
                output:
                tuple val(meta), path("*.bam"), emit: bam_ch

                script:
                "bwa mem ref reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)
        def yaml = MetaYmlGenerator.render(meta)
        def parsed = new Yaml().load(yaml) as Map

        then:
        // Named tuple: {bam_ch: [{meta: ...}, {'*.bam': ...}]}
        def output = parsed['output'] as List
        output.size() == 1
        output[0] instanceof Map
        def tupleItems = output[0]['bam_ch'] as List
        tupleItems.size() == 2
        tupleItems[0] == [meta: [type: 'map', description: 'TODO: Add description']]
        tupleItems[1] == ['*.bam': [type: 'list', description: 'TODO: Add description']]
    }

    def 'should render anonymous V1 tuple as nested YAML list'() {
        given:
        def mainNf = tempDir.resolve('main.nf')
        mainNf.text = '''\
            process ALIGN {
                input:
                tuple val(meta), path(reads)

                script:
                "bwa mem ref $reads"
            }
            '''.stripIndent()

        when:
        def meta = MetaYmlGenerator.extract(mainNf)
        def yaml = MetaYmlGenerator.render(meta)
        def parsed = new Yaml().load(yaml) as Map

        then:
        // Anonymous tuple (no emit name): rendered as nested list
        def input = parsed['input'] as List
        input.size() == 1
        input[0] instanceof List
        def tupleItems = input[0] as List
        tupleItems.size() == 2
        tupleItems[0] == [meta: [type: 'map', description: 'TODO: Add description']]
        tupleItems[1] == [reads: [type: 'file', description: 'TODO: Add description']]
    }

    // =========================================================================
    // Render tests
    // =========================================================================

    def 'should render YAML with comment header'() {
        given:
        def meta = new MetaYmlGenerator.ProcessMetadata(
            name: 'FASTQC',
            syntaxVersion: 'V1'
        )

        when:
        def yaml = MetaYmlGenerator.render(meta)

        then:
        yaml.startsWith('# This file was auto-generated by `nextflow module generate-meta`.')
        yaml.contains('# Review and complete all fields marked "TODO" before publishing.')
        yaml.contains("# Auto-detected types to verify:")
        yaml.contains("val: inferred as 'string'")
        yaml.contains("path: glob patterns")
    }

    def 'should omit input and output sections when empty'() {
        given:
        def meta = new MetaYmlGenerator.ProcessMetadata(
            name: 'SIMPLE',
            syntaxVersion: 'V1'
        )

        when:
        def yaml = MetaYmlGenerator.render(meta)

        then:
        !yaml.contains('input:')
        !yaml.contains('output:')
    }

    def 'should use RenderOptions values when provided'() {
        given:
        def meta = new MetaYmlGenerator.ProcessMetadata(name: 'MY_TOOL', syntaxVersion: 'V1')
        def opts = new MetaYmlGenerator.RenderOptions(
            name: 'nf-core/my-tool',
            version: '1.2.3',
            description: 'Does something useful',
            license: 'MIT',
            authors: ['Alice', 'Bob']
        )

        when:
        def yaml = MetaYmlGenerator.render(meta, opts)
        def parsed = new Yaml().load(yaml) as Map

        then:
        parsed['name'] == 'nf-core/my-tool'
        parsed['version'] == '1.2.3'
        parsed['description'] == 'Does something useful'
        parsed['license'] == 'MIT'
        parsed['authors'] == ['Alice', 'Bob']
    }

    def 'should fall back to TODO placeholders when RenderOptions fields are null'() {
        given:
        def meta = new MetaYmlGenerator.ProcessMetadata(name: 'MY_TOOL', syntaxVersion: 'V1')

        when:
        def yaml = MetaYmlGenerator.render(meta)
        def parsed = new Yaml().load(yaml) as Map

        then:
        (parsed['name'] as String).contains('my_tool')
        parsed['version'] == 'TODO: Add version'
        parsed['description'] == 'TODO: Add module description'
        (parsed['license'] as String).startsWith('TODO:')
        parsed['authors'] == ['TODO: Add author']
    }
}
