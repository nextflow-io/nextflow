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

package nextflow.script.parser.v2

import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptMeta
import test.Dsl2Spec

/**
 * Verify the compile-time collection of the {@code task} directives referenced by a
 * process definition.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDirectiveRefsTest extends Dsl2Spec {

    private Set<String> directiveRefs(String text, String processName='foo') {
        def parser = new ScriptLoaderV2(new Session())
        parser.setModule(true)
        parser.parse(text)
        parser.runScript()
        final process = ScriptMeta.get(parser.getScript()).getProcess(processName)
        return process.@taskBody.getDirectiveRefs()
    }

    def 'should collect the memory reference in the process script' () {
        expect:
        directiveRefs('''
            process foo {
                script:
                """
                java -Xmx${task.memory.toGiga()}g -jar app.jar
                """
            }
            ''') == ['memory'] as Set
    }

    def 'should collect the memory reference in a process directive' () {
        expect:
        directiveRefs('''
            process foo {
                beforeScript "echo ${task.memory}"

                script:
                """
                echo hello
                """
            }
            ''') == ['memory'] as Set
    }

    def 'should collect the memory reference in the process stub' () {
        expect:
        directiveRefs('''
            process foo {
                script:
                """
                echo hello
                """

                stub:
                """
                echo ${task.memory}
                """
            }
            ''') == ['memory'] as Set
    }

    def 'should collect the memory reference in the process ext directive' () {
        expect:
        directiveRefs('''
            process foo {
                ext args: "-Xmx${task.memory.toGiga()}g"

                script:
                """
                java ${task.ext.args} -jar app.jar
                """
            }
            ''') == ['memory', 'ext'] as Set
    }

    def 'should collect the memory reference behind a dynamic directive closure' () {
        expect:
        directiveRefs('''
            process foo {
                containerOptions { "--memory ${task.memory.toMega()}m" }

                script:
                """
                echo hello
                """
            }
            ''') == ['memory'] as Set
    }

    def 'should collect any referenced directive, not only the memory' () {
        expect:
        directiveRefs('''
            process foo {
                cpus 4

                script:
                """
                samtools sort -@ ${task.cpus} in.bam
                """
            }
            ''') == ['cpus'] as Set
    }

    def 'should not collect task properties that are not directives' () {
        expect:
        directiveRefs('''
            process foo {
                script:
                """
                echo "attempt ${task.attempt} of ${task.name} -- ${task.exitStatus}"
                """
            }
            ''') == [] as Set
    }

    def 'should report the referenced directive through the task run' () {
        given:
        def parser = new ScriptLoaderV2(new Session())
        parser.setModule(true)
        parser.parse('''
            process foo {
                script:
                """
                java -Xmx${task.memory.toGiga()}g -jar app.jar
                """
            }
            ''')
        parser.runScript()
        and:
        def processor = Mock(TaskProcessor) {
            getTaskBody() >> ScriptMeta.get(parser.getScript()).getProcess('foo').@taskBody
        }
        def task = new TaskRun(processor: processor)

        expect:
        task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }
}
