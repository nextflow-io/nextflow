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

package nextflow.script.control

import spock.lang.Specification
import test.TestUtils

import static test.TestUtils.deleteDir
import static test.TestUtils.tempDir
import static test.TestUtils.tempFile

/**
 * Type checking of calls to an included pipeline.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 */
class PipelineTypeCheckingTest extends Specification {

    static final String PIPELINE = '''\
        nextflow.enable.types = true

        params {
            input: Channel<String>
            fasta: Path
            aligner: String = 'star'
        }

        workflow {
            main:
            ch_bams = params.input.map { s -> file(s) }

            publish:
            bams = ch_bams
            multiqc = params.fasta
        }

        output {
            bams: Channel<Path> {}
            multiqc: Path {}
        }
        '''

    List<String> check(String main, Map<String,String> modules = [:]) {
        def root = tempDir()
        try {
            def mainFile = tempFile(root, 'main.nf', "nextflow.enable.types = true\n" + main.stripIndent())
            def pipelineFile = tempFile(root, 'rnaseq.nf', PIPELINE)
            def moduleFiles = modules.collect { name, text -> tempFile(root, name, "nextflow.enable.types = true\n" + text.stripIndent()) }
            def parser = new ScriptParser(root)
            return TestUtils.check(parser, [mainFile, pipelineFile] + moduleFiles)
                .findAll { e -> e.getSourceLocator().endsWith('main.nf') }
                .collect { e -> e.getOriginalMessage() }
        }
        finally {
            deleteDir(root)
        }
    }

    def 'should accept a channel argument for a channel param' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                RNASEQ( input: channel.of('a'), fasta: file('genome.fa') )
            }
            ''') == []
    }

    def 'should reject a dataflow argument for a param that is not a Channel or Value' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                RNASEQ( input: channel.of('a'), fasta: channel.value(file('genome.fa')) )
            }
            ''') == [ 'Param `fasta` expects a Path but received a Value<Path>' ]
    }

    def 'should reject a value argument for a channel param' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                RNASEQ( input: channel.value('a'), fasta: file('genome.fa') )
            }
            ''') == [ 'Param `input` expects a Channel<String> but received a Value<String>' ]
    }

    def 'should check the fields of a record argument' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                RNASEQ( record(input: channel.of('a'), aligner: 42) )
            }
            ''') == [
                'Pipeline `RNASEQ` requires the following params: fasta',
                'Param `aligner` expects a String but received a Integer'
            ]
    }

    def 'should check an included params record with overrides' () {
        expect:
        check('''\
            include { params as RnaseqParams ; workflow as RNASEQ } from './rnaseq.nf'

            params {
                rnaseq: RnaseqParams
            }

            workflow {
                RNASEQ( params.rnaseq + record(aligner: 42) )
            }
            ''') == [ 'Param `aligner` expects a String but received a Integer' ]
    }

    def 'should treat the fields of an included params record as nullable' () {
        expect:
        check('''\
            include { params as RnaseqParams } from './rnaseq.nf'

            workflow {
                main:
                SUMMARY( record(aligner: 'hisat2') )
                def p = record(aligner: 'hisat2') as RnaseqParams
            }

            workflow SUMMARY {
                take:
                p: RnaseqParams

                main:
                p.aligner
            }
            ''') == []
    }

    def 'should not require a partial record type param' () {
        expect:
        check('''\
            include { workflow as META } from './meta.nf'

            workflow {
                META( label: 'x' )
            }
            ''', [
            'meta.nf': '''\
            include { params as RnaseqParams ; workflow as RNASEQ } from './rnaseq.nf'

            params {
                rnaseq: RnaseqParams
                options: Options
                sample: Sample
                label: String
            }

            record Options {
                verbose: Boolean?
            }

            record Sample {
                id: String
            }

            workflow {
                RNASEQ( params.rnaseq )
            }
            '''
        ]) == [ 'Pipeline `META` requires the following params: sample' ]
    }

    def 'should report an invalid pipeline call' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                ''' + CALL + '''
            }
            ''') == [ ERROR ]

        where:
        CALL                                                                    | ERROR
        "RNASEQ( input: channel.of('a'), fasta: file('x'), foo: 1 )"            | 'Param `foo` is not defined by pipeline `RNASEQ`'
        "RNASEQ( input: channel.of('a') )"                                      | 'Pipeline `RNASEQ` requires the following params: fasta'
        "RNASEQ()"                                                              | 'Pipeline `RNASEQ` requires the following params: input, fasta'
        "RNASEQ( channel.of('a'), file('x') )"                                  | 'Pipeline `RNASEQ` should be called with named arguments, one for each of its params'
        "RNASEQ( 'a' )"                                                         | 'Pipeline `RNASEQ` should be called with named arguments or a record, but received a String'
    }

    def 'should report an included pipeline that is not aliased' () {
        expect:
        check('''\
            include { workflow ; params as RnaseqParams } from './rnaseq.nf'
            include { output } from './rnaseq.nf'
            ''') == [
                'An included pipeline must be aliased, e.g. `workflow as MY_PIPELINE`',
                'An included pipeline must be aliased, e.g. `output as MY_PIPELINE`'
            ]
    }

    def 'should return nothing from a pipeline without an output block' () {
        expect:
        check('''\
            include { workflow as HELLO } from './hello.nf'

            workflow {
                HELLO().foo
            }
            ''', [
            'hello.nf': '''\
            workflow {
                println 'Hello'
            }
            '''
        ]) == [ 'Unrecognized property `foo` for type void' ]
    }

    def 'should not wrap channel outputs' () {
        expect:
        check('''\
            include { workflow as RNASEQ } from './rnaseq.nf'

            workflow {
                r = RNASEQ( input: channel.of('a'), fasta: file('genome.fa') )
                BAMS( r.bams )
            }

            workflow BAMS {
                take:
                bams: Channel<Path>

                main:
                bams.view()
            }
            ''') == []
    }

    def 'should wrap non-channel outputs in a Value' () {
        expect:
        check('''\
            include { workflow as RNASEQ ; output as RnaseqOutput } from './rnaseq.nf'

            workflow {
                r = RNASEQ( input: channel.of('a'), fasta: file('genome.fa') )
                r.multiqc.view()
                SUMMARY( r )
            }

            workflow SUMMARY {
                take:
                r: RnaseqOutput

                main:
                r.multiqc.map { p -> p.name }.view()
            }
            ''') == []
    }

}
