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

import java.nio.file.Files

import nextflow.exception.ScriptCompilationException
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
}
