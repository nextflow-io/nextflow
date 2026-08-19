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
import nextflow.processor.TaskConfig
import nextflow.processor.TaskContext
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptMeta
import test.Dsl2Spec

/**
 * Verify that the {@code task} directives referenced by a process script are reported
 * through the variable references already collected for the task hash.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
class ProcessDirectiveRefsTest extends Dsl2Spec {

    private Set<String> valNames(String text, String processName='foo') {
        def parser = new ScriptLoaderV2(new Session())
        parser.setModule(true)
        parser.parse(text)
        parser.runScript()
        final process = ScriptMeta.get(parser.getScript()).getProcess(processName)
        return process.@taskBody.getValNames()
    }

    def 'should collect any referenced task property' () {
        expect:
        valNames('''
            process foo {
                script:
                """
                samtools sort -@ ${task.cpus} -m ${task.memory.toMega()}M in.bam
                """
            }
            ''') == ['task.cpus', 'task.memory'] as Set
    }

    def 'should collect the task ext and params references' () {
        expect:
        valNames('''
            process foo {
                script:
                """
                echo ${task.ext.args} ${params.foo}
                """
            }
            ''') == ['task.ext.args', 'params.foo'] as Set
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
            getDeclaredNames() >> []
        }
        def task = new TaskRun(processor: processor, context: new TaskContext(processor))

        expect:
        task.isDirectiveReferenced('memory')
        !task.isDirectiveReferenced('cpus')
    }

    def 'should not add the task references to the task hash' () {
        given:
        // `task` is put in the task context by TaskConfig#setContext, therefore every
        // `task.*` reference is skipped by #getGlobalVars as a task-local variable. This is
        // what allows the directive references to travel on `valRefs` without invalidating
        // the resume cache of every pipeline -- `task.ext.*` is the documented exception,
        // re-added to the hash by TaskHasher#getTaskExtensionDirectiveVars.
        def context = new TaskContext(holder: [:])
        def config = new TaskConfig(memory: '8 GB').setContext(context)
        and:
        def task = Spy(TaskRun)
        task.processor = Mock(TaskProcessor) { getName() >> 'foo' }
        task.context = context
        task.config = config
        task.getVariableNames() >> (['task.memory', 'task.ext.args', 'x'] as Set)

        when:
        def vars = task.getGlobalVars(new Binding(x: 1))

        then:
        context.isLocalVar('task')
        and:
        vars == [x: 1]
    }
}
