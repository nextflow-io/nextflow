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

import test.Dsl2Spec

import static test.ScriptHelper.*

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class DataflowTypesTest extends Dsl2Spec {

    def 'should allow legacy workflow to call legacy process/workflow' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            include { foo ; bar } from './module.nf'

            workflow {
                ch = channel.of(1, 2, 3)
                bar(foo(ch))
            }
            '''

        folder.resolve('module.nf').text = '''
            process foo {
                input:
                val x

                output:
                val y

                exec:
                y = x * 2
            }

            workflow bar {
                take:
                ys

                emit:
                ys.map { y -> y % 3 }
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should allow legacy workflow to call typed process/workflow' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            include { foo ; bar } from './module.nf'

            workflow {
                ch = channel.of(1, 2, 3)
                bar(foo(ch))
            }
            '''

        folder.resolve('module.nf').text = '''
            nextflow.enable.types = true

            process foo {
                input:
                x: Integer

                output:
                y: Integer

                exec:
                y = x * 2
            }

            workflow bar {
                take:
                ys: Channel<Integer>

                emit:
                ys.map { y -> y * 2 }
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should allow typed workflow to call legacy process/workflow' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { foo ; bar } from './module.nf'

            workflow {
                ch = channel.of(1, 2, 3)
                bar(foo(ch))
            }
            '''

        folder.resolve('module.nf').text = '''
            process foo {
                input:
                val x

                output:
                val y

                exec:
                y = x * 2
            }

            workflow bar {
                take:
                ys

                emit:
                ys.map { y -> y % 3 }
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should allow typed workflow to call typed process/workflow' () {
        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { foo ; bar } from './module.nf'

            workflow {
                ch = channel.of(1, 2, 3)
                bar(foo(ch))
            }
            '''

        folder.resolve('module.nf').text = '''
            nextflow.enable.types = true

            process foo {
                input:
                x: Integer

                output:
                y: Integer

                exec:
                y = x * 2
            }

            workflow bar {
                take:
                ys: Channel<Integer>

                emit:
                ys.map { y -> y * 2 }
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
