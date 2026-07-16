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

import java.lang.reflect.ParameterizedType
import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.ScriptCompilationException
import nextflow.util.RecordMap
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptTypesTest extends Dsl2Spec {

    def 'should support enum types' () {

        when:
        def result = runScript(
            '''\
            enum Day {
                MONDAY,
                TUESDAY,
                WEDNESDAY,
                THURSDAY,
                FRIDAY,
                SATURDAY,
                SUNDAY
            }

            workflow {
                Day.TUESDAY
            }
            '''
        )
        then:
        result.toString() == 'TUESDAY'
    }

    def 'should support record types' () {

        when:
        def result = runScript(
            '''\
            workflow {
                s = record(id: '1', fastq: file('1.fastq'))
                printSample(s)
            }

            def printSample(s: Sample) {
                return "Sample(id: ${s.id}, fastq: ${s.fastq.name})"
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''
        )
        then:
        result == 'Sample(id: 1, fastq: 1.fastq)'
    }

    def 'should cast value to record type' () {

        when:
        def result = runScript(
            '''\
            workflow {
                record(id: '1', fastq: '1.fastq') as Sample
            }

            record Sample {
                id: String
                fastq: String
            }
            '''
        )
        then:
        result instanceof RecordMap
        result.id == '1'
        result.fastq == '1.fastq'
    }

    def 'should cast value to parameterized type'() {
        when:
        def samples = runScript(
            '''\
            workflow {
                [
                  [id: '1', fastq: '1.fastq'],
                  [id: '2', fastq: '2.fastq'],
                  [id: '3', fastq: '3.fastq']
                ] as List<Sample>
            }

            record Sample {
                id: String
                fastq: String
            }
            '''
        )
        then:
        samples instanceof List
        samples.size() == 3
        samples[0] instanceof RecordMap
        samples[0].id == '1'
        samples[0].fastq == '1.fastq'
        samples[1] instanceof RecordMap
        samples[1].id == '2'
        samples[1].fastq == '2.fastq'
        samples[2] instanceof RecordMap
        samples[2].id == '3'
        samples[2].fastq == '3.fastq'
    }

    def 'should strip unsupported type annotations' () {

        when:
        runScript(
            '''\
            // strip Tuple type from cast expression
            ['1', '1.fastq', '2.fastq'] as Tuple<String,String,String>

            // strip type annotation in variable declaration
            def ch: Channel = channel.empty()
            '''
        )

        then:
        noExceptionThrown()
    }

    def 'should replace instanceof on record type with runtime check' () {

        when:
        def result = runScript(
            '''\
            workflow {
                [
                    test1: record(id: '1', fastq: file('1.fastq')) instanceof Sample,
                    test2: record(id: '2') instanceof Sample,
                    test3: 42 instanceof Integer,
                    test4: '42' instanceof Integer
                ]
            }

            record Sample {
               	id: String
               	fastq: Path
            }
            '''
        )
        then:
        result.test1 == true
        result.test2 == false
        result.test3 == true
        result.test4 == false
    }

    def 'should report error for invalid record call' () {

        when:
        runScript(
            '''\
            workflow {
                Sample(id: '1', fastq: file('1.fastq'))
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''
        )
        then:
        def e = thrown(ScriptCompilationException)
        e.cause.message.contains '`Sample` is not defined'
    }

    def 'should allow record types to be included across modules' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            include { Sample } from './module.nf'

            workflow {
                s = record(id: '1', fastq: file('1.fastq'))
                printSample(s)
            }

            def printSample(s: Sample) {
                println s
            }
            '''

        folder.resolve('module.nf').text = '''
            record Sample {
                id: String
                fastq: Path
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should allow enum types to be included across modules' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            include { Color } from './module.nf'

            workflow {
                println Color.RED
            }
            '''

        folder.resolve('module.nf').text = '''
            enum Color {
                RED,
                GREEN,
                BLUE
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should allow record types with same name in different modules' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            include { printSample } from './module.nf'

            workflow {
                s = record(id: '1', fastq: file('1.fastq'))
                printSample(s)
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''

        folder.resolve('module.nf').text = '''
            def printSample(s: Sample) {
                println s
            }

            record Sample {
                id: String
                fastq: Path
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should stage path fields of an included record type' () {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('types.nf').text = '''\
            record Pair {
                id: String
                r1: Path
                r2: Path
            }
            '''.stripIndent()

        def file = folder.resolve('main.nf')
        file.text = '''\
            nextflow.enable.types = true

            include { Pair } from './types.nf'

            process CONSUME {
                input:
                rec: Pair

                script:
                """
                cat ${rec.r1} ${rec.r2}
                """
            }

            workflow {
            }
            '''.stripIndent()

        when:
        runScript(file)
        and:
        def process = ScriptMeta.get(ScriptMeta.getScriptByPath(file)).getProcess('CONSUME')
        def inputs = process.getProcessConfig().getInputs()

        then:
        // the `r1` and `r2` record fields should be staged as input files
        inputs.getFiles().size() == 2

        cleanup:
        folder.deleteDir()
    }

    def 'should expose type annotations via reflection'() {
        when:
        def script = loadScript(
            '''\
            record Sample {
                id: String
                reads: List<Path>
            }
            ''',
            module: true
        )
        def meta = ScriptMeta.get(script)
        def typeDef = meta.getComponent('Sample') as TypeDef
        def type = typeDef.getTarget()
        then:
        type.getField('id').getType() == String
        type.getField('id').getGenericType() instanceof Class
        type.getField('id').getGenericType() == String
        type.getField('reads').getType() == List
        type.getField('reads').getGenericType() instanceof ParameterizedType
        type.getField('reads').getGenericType().getRawType() == List
        type.getField('reads').getGenericType().getActualTypeArguments()[0] instanceof Class
        type.getField('reads').getGenericType().getActualTypeArguments()[0] == Path

        when:
        script = loadScript(
            '''\
            def greet(message: String, names: List<String>) {
            }
            ''',
            module: true
        )
        def method = script.getClass().getDeclaredMethods().find { m -> m.name == 'greet' }
        def params = method.getParameters()
        then:
        params[0].getType() == String
        params[0].getParameterizedType() instanceof Class
        params[0].getParameterizedType() == String
        params[1].getType() == List
        params[1].getParameterizedType() instanceof ParameterizedType
        params[1].getParameterizedType().getRawType() == List
        params[1].getParameterizedType().getActualTypeArguments()[0] instanceof Class
        params[1].getParameterizedType().getActualTypeArguments()[0] == String
    }
}
