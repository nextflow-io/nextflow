/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.Path

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractGridExecutorTest extends Specification {

    def 'should remove invalid chars from name' () {

        given:
        def task = new TaskRun(name: 'task 90 = (foo:bar/baz)')
        def exec = [:] as AbstractGridExecutor

        expect:
        exec.getJobNameFor(task) == 'nf-task_90___(foo_bar_baz)'

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

    def 'should cut long job name' () {

        given:
        def exec = [:] as AbstractGridExecutor
        def task = Mock(TaskRun)

        when:
        def name = exec.getJobNameFor(task)
        then:
        1 * task.getName() >> { 'abcd' * 100 }
        and:
        name.size() == 256
        name == ( 'nf-' + 'abcd' * 100 ).substring(0,256)
    }

    def 'should sanitize job name' () {
        given:
        def LONG = 'abcd' * 100
        def exec = [:] as AbstractGridExecutor
        
        expect:
        exec.sanitizeJobName('foo') == 'foo'
        exec.sanitizeJobName(LONG) == LONG.substring(0,256)
    }

    def 'should add change dir variable' () {
        given:
        def work = Path.of('/some/dir')
        def exec = Spy(AbstractGridExecutor)
        def task = Mock(TaskRun) { getWorkDir() >> work}
        when:
        def result = exec.getHeaderScript(task)
        then:
        1 * exec.getHeaders(task) >> '#$ one\n#$ two\n'
        result == '''\
                #$ one
                #$ two
                NXF_CHDIR=/some/dir
                '''.stripIndent()
    }
    
    def 'should fetch queue status'() {
        given:
        def STATUS = ['123': AbstractGridExecutor.QueueStatus.RUNNING]
        def NAME = 'TheExecutorName'
        and:
        def session = Mock(Session) { getConfig()>>[:] }
        and:
        def exec = Spy(AbstractGridExecutor)
        exec.session = session
        exec.@queueInterval = Duration.of('1m')
        exec.name = NAME


        when:
        def result = exec.getQueueStatus('foo')
        then:
        1 * session.getExecConfigProp(NAME,'queueGlobalStatus',false)>>false
        1 * exec.getQueueStatus0('foo') >> STATUS
        and:
        result == STATUS


        when:
        result = exec.getQueueStatus('foo')
        then:
        1 * session.getExecConfigProp(NAME,'queueGlobalStatus',false)>>true
        1 * exec.getQueueStatus0(null) >> STATUS
        and:
        result == STATUS
    }

    def 'should add cluster options' () {
        given:
        def exec = Spy(AbstractGridExecutor)

        when:
        def result = []
        exec.addClusterOptionsDirective(new TaskConfig(clusterOptions: OPTS), result)
        then:
        result == EXPECTED

        where:
        OPTS                    | EXPECTED
        null                    | []
        '-foo 1'                | ['-foo 1', '']
        '-foo 1 --bar 2'        | ['-foo 1 --bar 2', '']
        ['-foo 1','--bar 2']    | ['-foo 1', '', '--bar 2', '']
    }

}
