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
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
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
        meta.outputs[0].resolvedType == 'file'
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
        meta.outputs[0].subEntries[1].resolvedType == 'file'
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

        then:
        yaml.startsWith('# This file was auto-generated')
        yaml.contains('name:')
        yaml.contains('tools:')
        yaml.contains('input:')
        yaml.contains('output:')
        yaml.contains('meta:')
        yaml.contains('reads:')
        yaml.contains("type: map")
        yaml.contains("type: file")
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
        meta.inputs[2].resolvedType == 'int'
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
        meta.inputs.size() == 3
        meta.inputs[0].isTuple()
        meta.inputs[0].subEntries.size() == 2
        meta.inputs[0].subEntries[0].identifier == 'sample_id'
        meta.inputs[0].subEntries[0].resolvedType == 'String'
        meta.inputs[0].subEntries[1].identifier == 'reads'
        meta.inputs[0].subEntries[1].resolvedType == 'file'
        meta.inputs[1].identifier == 'index'
        meta.inputs[1].resolvedType == 'file'
        meta.inputs[2].identifier == 'id'
        meta.inputs[2].resolvedType == 'string'

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
        meta.outputs[1].resolvedType == 'int'
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
        yaml.contains('# Fields marked "TODO" require manual completion before publishing.')
        yaml.contains('# Required fields: name, version, description, license, authors')
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

    def 'should include tools placeholder always'() {
        given:
        def meta = new MetaYmlGenerator.ProcessMetadata(
            name: 'MY_TOOL',
            syntaxVersion: 'V1'
        )

        when:
        def yaml = MetaYmlGenerator.render(meta)

        then:
        yaml.contains('tools:')
        yaml.contains('TODO-tool-name')
        yaml.contains('TODO: Add tool description')
    }
}
