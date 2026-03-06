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
 * Tests for ModuleSpec
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleSpecTest extends Specification {

    @TempDir
    Path tempDir

    def 'should load valid manifest' () {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
        metaYaml.text = '''
name: nf-core/fastqc
version: 1.0.0
description: FastQC quality control
authors:
  - John Doe
license: MIT
keywords:
  - quality-control
  - fastq
requires:
  nextflow: ">=24.04.0"
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.name == 'nf-core/fastqc'
        manifest.version == '1.0.0'
        manifest.description == 'FastQC quality control'
        manifest.authors == ['John Doe']
        manifest.license == 'MIT'
        manifest.keywords == ['quality-control', 'fastq']
        manifest.requires == ['nextflow': '>=24.04.0']
    }

    def 'should fail to load non-existent manifest' () {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')

        when:
        ModuleSpec.load(metaYaml)

        then:
        thrown(AbortOperationException)
    }

    def 'should load manifest with simple inputs and outputs'() {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
        metaYaml.text = '''\
name: nf-core/fastqc
version: 1.0.0
description: FastQC quality control
license: MIT
input:
  - sample_id:
      type: string
      description: Sample identifier
  - reads:
      type: file
      description: Input reads file
output:
  - html:
      type: list
      description: FastQC HTML reports
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.inputs.size() == 2
        manifest.inputs[0].identifier == 'sample_id'
        manifest.inputs[0].type == 'string'
        manifest.inputs[0].description == 'Sample identifier'
        manifest.inputs[1].identifier == 'reads'
        manifest.inputs[1].type == 'file'
        manifest.inputs[1].description == 'Input reads file'

        and:
        manifest.outputs.size() == 1
        manifest.outputs[0].identifier == 'html'
        manifest.outputs[0].type == 'list'
        manifest.outputs[0].description == 'FastQC HTML reports'
    }

    def 'should load manifest with anonymous tuple input'() {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
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
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.inputs.size() == 1
        manifest.inputs[0].isTuple()
        manifest.inputs[0].identifier == 'tuple'
        manifest.inputs[0].subEntries.size() == 2
        manifest.inputs[0].subEntries[0].identifier == 'meta'
        manifest.inputs[0].subEntries[0].type == 'map'
        manifest.inputs[0].subEntries[1].identifier == 'reads'
        manifest.inputs[0].subEntries[1].type == 'file'
    }

    def 'should load manifest with named tuple output'() {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
        metaYaml.text = '''\
name: nf-core/align
version: 1.0.0
description: Align reads
license: MIT
output:
  - bam_ch:
    - meta:
        type: map
        description: Sample metadata
    - bam:
        type: file
        description: Output BAM file
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.outputs.size() == 1
        manifest.outputs[0].isTuple()
        manifest.outputs[0].identifier == 'bam_ch'
        manifest.outputs[0].subEntries.size() == 2
        manifest.outputs[0].subEntries[0].identifier == 'meta'
        manifest.outputs[0].subEntries[0].type == 'map'
        manifest.outputs[0].subEntries[1].identifier == 'bam'
        manifest.outputs[0].subEntries[1].type == 'file'
    }

    def 'should load manifest with missing type using TODO placeholder'() {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
        metaYaml.text = '''\
name: nf-core/fastqc
version: 1.0.0
description: FastQC quality control
license: MIT
input:
  - reads:
      description: Input reads
output:
  - report:
      type: file
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.inputs[0].identifier == 'reads'
        manifest.inputs[0].type == ModuleSpec.TYPE_NOT_FOUND_TODO
        manifest.outputs[0].identifier == 'report'
        manifest.outputs[0].description == ModuleSpec.DESCRIPTION_TODO
    }

    def 'should load manifest with no inputs or outputs'() {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
        metaYaml.text = '''\
name: nf-core/fastqc
version: 1.0.0
description: FastQC quality control
license: MIT
'''

        when:
        def manifest = ModuleSpec.load(metaYaml)

        then:
        manifest.inputs.isEmpty()
        manifest.outputs.isEmpty()
    }

    def 'should validate complete manifest' () {
        given:
        def manifest = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: '1.0.0',
            description: 'FastQC quality control',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

        then:
        errors.isEmpty()
        manifest.isValid()
    }

    def 'should detect missing required fields' () {
        given:
        def manifest = new ModuleSpec(
            name: 'nf-core/fastqc'
            // missing version, description, license
        )

        when:
        def errors = manifest.validate()

        then:
        errors.size() == 3
        errors.any { it.contains('version') }
        errors.any { it.contains('description') }
        errors.any { it.contains('license') }
        !manifest.isValid()
    }

    def 'should validate version format' () {
        given:
        def manifest = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: version,
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

        then:
        errors.isEmpty() == valid

        where:
        version         | valid
        '1.0.0'         | true
        '1.0.0-alpha'   | true
        '1.0.0-beta.1'  | true
        '1.0'           | false
        'v1.0.0'        | false
        '1.0.0.0'       | false
    }

    def 'should validate module name format' () {
        given:
        def manifest = new ModuleSpec(
            name: name,
            version: '1.0.0',
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

        then:
        errors.isEmpty() == valid

        where:
        name                        | valid
        'nf-core/fastqc'            | true
        'myorg/my-module'           | true
        'org_1/tool_2'              | true
        'nf-core/gfatools/gfa2fa'   | true   // nested module path
        'myorg/tools/sub/module'    | true   // deeply nested
        'org.name/tool/sub'         | true   // dot in scope
        'fastqc'                    | false
        '@nf-core/fastqc'           | false
        'nf-core/fast qc'           | false
        'nf-core/'                  | false  // trailing slash
        '/nf-core/fastqc'           | false  // leading slash
    }

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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.name == 'my-scope/fastqc'
        spec.inputs.size() == 2
        spec.inputs[0].identifier == 'sample_id'
        spec.inputs[0].type == 'string'
        spec.inputs[1].identifier == 'reads'
        spec.inputs[1].type == 'file'
        spec.outputs.size() == 1
        spec.outputs[0].identifier == '*.html'
        spec.outputs[0].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 1
        spec.inputs[0].isTuple()
        spec.inputs[0].subEntries.size() == 2
        spec.inputs[0].subEntries[0].identifier == 'meta'
        spec.inputs[0].subEntries[0].type == 'map'
        spec.inputs[0].subEntries[1].identifier == 'reads'
        spec.inputs[0].subEntries[1].type == 'file'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].isTuple()
        spec.outputs[0].subEntries.size() == 2
        spec.outputs[0].subEntries[0].identifier == 'meta'
        spec.outputs[0].subEntries[0].type == 'map'
        spec.outputs[0].subEntries[1].identifier == '*.bam'
        spec.outputs[0].subEntries[1].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.isEmpty()
        spec.outputs.size() == 1
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))
        def yaml = spec.render()

        then:
        spec.name == 'my-scope/samtools_sort'
        yaml.contains('my-scope/samtools_sort')
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))
        def yaml = spec.render()
        def parsed = new Yaml().load(yaml) as Map

        then:
        // Top-level structure
        yaml.startsWith('# This file was auto-generated')
        parsed['name'] == 'my-scope/fastqc'
        // null fields are omitted from the YAML map; TODO list is in the comment header
        !parsed.containsKey('version')
        yaml.contains('# TODO List:')
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.name == 'my-scope/first'
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
        ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 3
        spec.inputs[0].identifier == 'sample_id'
        spec.inputs[0].type == 'string'
        spec.inputs[1].identifier == 'reads'
        spec.inputs[1].type == 'file'
        spec.inputs[2].identifier == 'index'
        spec.inputs[2].type == 'integer'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bam'
        spec.outputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        // 4 inputs: tuple(sample_id, reads), index, id, anon_tuple
        spec.inputs.size() == 4
        spec.inputs[0].isTuple()
        spec.inputs[0].subEntries.size() == 2
        spec.inputs[0].subEntries[0].identifier == 'sample_id'
        spec.inputs[0].subEntries[0].type == 'string'
        spec.inputs[0].subEntries[1].identifier == 'reads'
        spec.inputs[0].subEntries[1].type == 'file'
        spec.inputs[1].identifier == 'index'
        spec.inputs[1].type == 'file'
        spec.inputs[2].identifier == 'id'
        spec.inputs[2].type == 'string'
        !spec.inputs[3].isTuple()
        spec.inputs[3].identifier == 'anon_tuple'
        spec.inputs[3].type == 'tuple'

        and:
        // 4 outputs inferred via symbol map
        spec.outputs.size() == 4
        // bam = file("out.bam") → RHS is file() call → file
        spec.outputs[0].identifier == 'bam'
        spec.outputs[0].type == 'file'
        // out_tuple = tuple(sample_id, file("out2.bam")) → tuple with 2 sub-entries
        spec.outputs[1].identifier == 'out_tuple'
        spec.outputs[1].isTuple()
        spec.outputs[1].subEntries.size() == 2
        spec.outputs[1].subEntries[0].identifier == 'sample_id'
        spec.outputs[1].subEntries[0].type == 'string'   // from input symbol map
        spec.outputs[1].subEntries[1].identifier == 'out2.bam'  // from file("out2.bam") arg
        spec.outputs[1].subEntries[1].type == 'file'
        // res = test → exec: test = "hola" → symbol map: test → string
        spec.outputs[2].identifier == 'res'
        spec.outputs[2].type == 'string'
        // res2 = id → input id: String → symbol map: id → string
        spec.outputs[3].identifier == 'res2'
        spec.outputs[3].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs[0].identifier == 'typed_input'
        spec.inputs[0].type == 'file'
        // exec: result = test; test = 'done' → symbol map has result → string
        // so exec block assignment overrides the tuple() in the output block
        spec.outputs.size() == 2
        spec.outputs[0].identifier == 'result'
        spec.outputs[0].type == 'string'
        // outIteration = iteration → declared type int → integer
        spec.outputs[1].identifier == 'outIteration'
        spec.outputs[1].type == 'integer'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'sorted.bam'
        spec.outputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].type == 'list'

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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 1
        spec.inputs[0].identifier == 'reads'
        spec.inputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 1
        spec.inputs[0].identifier == 'GENOME_DIR'
        spec.inputs[0].type == 'string'

        and:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'ALIGNED_READS'
        spec.outputs[0].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 1
        spec.inputs[0].identifier == 'stdin'
        spec.inputs[0].type == 'stdin'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'stdout'
        spec.outputs[0].type == 'stdout'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'count'
        spec.outputs[0].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs.size() == 1
        spec.inputs[0].identifier == 'reads'
        spec.inputs[0].type == 'file'

        and:
        spec.outputs.size() == 2
        spec.outputs[0].identifier == 'result.txt'
        spec.outputs[0].type == 'file'
        spec.outputs[1].identifier == '*.log'
        spec.outputs[1].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.inputs[0].isTuple()
        spec.inputs[0].subEntries.size() == 3
        spec.inputs[0].subEntries[0].identifier == 'sample_id'
        spec.inputs[0].subEntries[0].type == 'string'
        spec.inputs[0].subEntries[1].identifier == 'reads'
        spec.inputs[0].subEntries[1].type == 'file'
        spec.inputs[0].subEntries[2].identifier == 'GENOME'
        spec.inputs[0].subEntries[2].type == 'string'

        and:
        spec.outputs[0].isTuple()
        spec.outputs[0].subEntries.size() == 3
        spec.outputs[0].subEntries[0].identifier == 'sample_id'
        spec.outputs[0].subEntries[0].type == 'string'
        spec.outputs[0].subEntries[1].identifier == '*.bam'
        spec.outputs[0].subEntries[1].type == 'list'
        spec.outputs[0].subEntries[2].identifier == 'ALIGNED_READS'
        spec.outputs[0].subEntries[2].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bam'
        spec.outputs[0].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bam'
        spec.outputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bam'
        spec.outputs[0].type == 'file'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bams'
        spec.outputs[0].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'result'
        spec.outputs[0].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'result'
        spec.outputs[0].type == 'string'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))

        then:
        spec.outputs.size() == 1
        spec.outputs[0].identifier == 'bam_ch'
        spec.outputs[0].isTuple()
        spec.outputs[0].subEntries[1].type == 'list'
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))
        def yaml = spec.render()
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
        def spec = ModuleSpec.extract(mainNf, new ModuleSpec.ExtractOptions(scope: 'my-scope'))
        def yaml = spec.render()
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
        def spec = new ModuleSpec(name: 'my-scope/fastqc')

        when:
        def yaml = spec.render()

        then:
        yaml.startsWith('# This file was auto-generated by `nextflow module spec`.')
        yaml.contains('# Review and complete all fields marked "TODO" before publishing.')
        yaml.contains("# Auto-detected types to verify:")
        yaml.contains("val: inferred as 'string'")
        yaml.contains("path: glob patterns")
    }

    def 'should render TODO list for missing required fields'() {
        given:
        def spec = new ModuleSpec(name: 'my-scope/fastqc')
        // version, description, license are all missing

        when:
        def yaml = spec.render()

        then:
        yaml.contains('# TODO List:')
        yaml.contains('Missing required field: version')
        yaml.contains('Missing required field: description')
        yaml.contains('Missing required field: license')
    }

    def 'should not render TODO list when all required fields are present'() {
        given:
        def spec = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: '1.0.0',
            description: 'Run FastQC',
            license: 'MIT'
        )

        when:
        def yaml = spec.render()

        then:
        !yaml.contains('# TODO List:')
    }

    def 'should omit null fields and include present fields in rendered YAML'() {
        given:
        def spec = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: '1.0.0',
            description: 'Run FastQC',
            license: 'MIT'
        )

        when:
        def parsed = new Yaml().load(spec.render()) as Map

        then:
        parsed['name'] == 'nf-core/fastqc'
        parsed['version'] == '1.0.0'
        parsed['description'] == 'Run FastQC'
        parsed['license'] == 'MIT'
    }

    def 'should omit input and output sections when empty'() {
        given:
        def spec = new ModuleSpec(name: 'my-scope/simple')

        when:
        def yaml = spec.render()

        then:
        !yaml.contains('input:')
        !yaml.contains('output:')
    }

}
