/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowVariable
import org.junit.Rule

import nextflow.Channel
import test.Dsl2Spec
import test.OutputCapture

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MultiMapOpTest extends Dsl2Spec {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should fork channel' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(0,1,2)
                .multiMap {
                        foo: it+1
                        bar: it*it+2
                        baz: 3
                    }
        ''')
        then:
        result.size() == 3
        and:
        result[0].val == 1
        result[0].val == 2
        result[0].val == 3
        result[0].val == Channel.STOP
        and:
        result[1].val == 2
        result[1].val == 3
        result[1].val == 6
        result[1].val == Channel.STOP
        and:
        result[2].val == 3
        result[2].val == 3
        result[2].val == 3
        result[2].val == Channel.STOP

    }

    def 'should fork channel with custom param' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(0,1,2)
                .multiMap { p ->
                        foo: p+1
                        bar: p*p+2
                        baz: p-1
                    }
        ''')
        then:
        result.size() == 3
        and:
        result[0].val == 1
        result[0].val == 2
        result[0].val == 3
        result[0].val == Channel.STOP
        and:
        result[1].val == 2
        result[1].val == 3
        result[1].val == 6
        result[1].val == Channel.STOP
        and:
        result[2].val == -1
        result[2].val ==  0
        result[2].val ==  1
        result[2].val == Channel.STOP

    }

    def 'should pass criteria as argument' () {
        when:
        dsl_eval('''   
            criteria = multiMapCriteria { 
                foo: it
                bar: it*it
            }

            ch1 = Channel.from(1,2,3).multiMap(criteria)  
            
            ch1.foo.view { "foo:$it" }
            ch1.bar.view { "bar:$it" }
        ''')

        def stdout = capture.toString()
        then:
        stdout.contains('foo:1')
        stdout.contains('foo:2')
        stdout.contains('foo:3')
        and:
        stdout.contains('bar:1')
        stdout.contains('bar:4')
        stdout.contains('bar:9')
    }

    def 'should fork channel value ch' () {

        when:
        def result = dsl_eval('''   
            Channel
                .value('hello')
                .multiMap { p ->
                        foo: p.toUpperCase()
                        bar: p.reverse()
                    }
        ''')
        then:
        result.size() == 2
        and:
        result[0] instanceof DataflowVariable
        result[0].val == 'HELLO'
        and:
        result[1] instanceof DataflowVariable
        result[1].val == 'olleh'

    }

}
