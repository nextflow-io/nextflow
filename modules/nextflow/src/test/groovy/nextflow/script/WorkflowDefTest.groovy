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

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.exception.ScriptCompilationException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Timeout(5)
class WorkflowDefTest extends Dsl2Spec {

    def 'should define named workflows' () {

        given:
        def SCRIPT = '''

            workflow alpha {
              print 'Hello world'
            }

            workflow bravo {
              take:
              foo
              bar

              main:
              print foo
              print bar

              emit:
              foo+bar
            }

            workflow delta {
                take:
                foo
                bar

                main:
                println foo+bar
            }

            workflow empty { }
        '''

        when:
        def script = loadScript(SCRIPT, module: true)
        def meta = ScriptMeta.get(script)

        then:
        meta.definitions.size() == 4
        meta.getWorkflow('alpha') .declaredInputs == []
        meta.getWorkflow('bravo') .declaredInputs == ['foo', 'bar']
        meta.getWorkflow('delta') .declaredInputs == ['foo','bar']
        meta.getWorkflow('empty') .declaredInputs == []
    }

    def 'should define entry workflow' () {

        def SCRIPT = '''

            workflow {
              print 1
              print 2
            }
        '''

        when:
        def script = loadScript(SCRIPT)
        def meta = ScriptMeta.get(script)
        then:
        meta.getWorkflow(null)

    }

    def 'should run workflow block' () {

        given:
        def SCRIPT = '''

            workflow alpha {
              take: foo
              main: bar = foo ; baz = foo
              emit: bar ; baz
            }

        '''

        when:
        def script = loadScript(SCRIPT, module: true)
        def workflow = ScriptMeta.get(script).getWorkflow('alpha')
        then:
        workflow.declaredInputs == ['foo']
        workflow.declaredOutputs == ['bar', 'baz']

    }

    def 'should report malformed workflow block' () {

        given:
        def SCRIPT = '''

            workflow alpha {
              take: foo
              main: println foo
              take: bar
            }

        '''

        when:
        loadScript(SCRIPT)
        then:
        def e = thrown(ScriptCompilationException)
        e.cause.message.contains('Invalid workflow definition')

    }

    def 'should not fail' () {
        given:
        // this test checks that the closure used to define the workflow
        // does NOT define an implicit `it` parameter that would clash
        // with the `it` used by the inner closure

        def SCRIPT = """

        workflow {
            channel.empty().map { id -> id +1 }
            channel.empty().map { it -> def id = it+1 }
        }
        """

        when:
        runScript(SCRIPT)

        then:
        noExceptionThrown()
    }

    def 'should validate collect output'() {
        given:
        def ch1 = new DataflowQueue(); ch1 << 'blah blah'
        def ch2 = new DataflowQueue(); ch2 << 'xxx'

        def binding = new WorkflowBinding(foo: 'Hello', bar: 'world', ch1: ch1, ch2: new ChannelOut([ch2]))
        def workflow = new WorkflowDef(binding: binding)

        when:
        def result = workflow.collectOutputs(['foo'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowVariable
        result.foo instanceof DataflowVariable
        result.foo.val == 'Hello'

        when:
        result = workflow.collectOutputs(['foo', 'bar'])
        then:
        result instanceof ChannelOut
        result.size()==2
        result[0] instanceof DataflowVariable
        result[1] instanceof DataflowVariable
        result.foo.val == 'Hello'
        result.bar.val == 'world'

        when:
        result = workflow.collectOutputs(['ch1'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowQueue
        result[0].val == 'blah blah'


        when:
        result = workflow.collectOutputs(['ch2'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowQueue
        result[0].val == 'xxx'

    }

    def 'should clone with a new name' () {
        given:
        def work = new WorkflowDef(name:'woo', body: new BodyDef({}, 'source'))

        when:
        def copy = work.cloneWithName('bar')
        then:
        copy.getName() == 'bar'
    }

}
