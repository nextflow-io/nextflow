/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import nextflow.Session
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractGridExecutorTest extends Specification {

    def 'should remove invalid chars from name' () {

        given:
        def task = new TaskRun(name: 'task 90 (foo:bar/baz)')
        def exec = [:] as AbstractGridExecutor

        expect:
        exec.getJobNameFor(task) == 'nf-task_90_(foo_bar_baz)'

    }

    def 'should return the kill list' () {

        given:
        def exec = [getKillCommand: { ['qdel'] }] as AbstractGridExecutor

        expect:
        exec.killTaskCommand('10') == ['qdel', '10']
        exec.killTaskCommand([11,12]) == ['qdel', '11', '12']
        exec.killTaskCommand([100,200,300]) == ['qdel', '100', '200', '300']

    }

    def 'should return a custom job name'() {

        given:
        def exec = [:] as AbstractGridExecutor
        exec.session = [:] as Session
        exec.session.config = [:]

        expect:
        exec.resolveCustomJobName(Mock(TaskRun)) == null

        when:
        exec.session = [:] as Session
        exec.session.config = [ executor: [jobName: { task.name.replace(' ','_') }  ] ]
        then:
        exec.resolveCustomJobName(new TaskRun(config: [name: 'hello world'])) == 'hello_world'

    }

    def 'should return job submit name' () {

        given:
        def exec = [:] as AbstractGridExecutor
        exec.session = [:] as Session
        exec.session.config = [:]

        final taskName = 'Hello world'
        final taskRun = new TaskRun(name: taskName, config: [name: taskName])

        expect:
        exec.getJobNameFor(taskRun) == 'nf-Hello_world'

        when:
        exec.session = [:] as Session
        exec.session.config = [ executor: [jobName: { task.name.replace(' ','_') }  ] ]
        then:
        exec.getJobNameFor(taskRun) == 'Hello_world'
    }
}
