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
import spock.lang.Ignore
import org.junit.Rule

import nextflow.Channel
import nextflow.exception.ScriptCompilationException
import test.Dsl2Spec
import test.OutputCapture
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BranchOpTest extends Dsl2Spec  {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should branch input values' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(0,1,2)
                .branch {
                        foo: it <1
                        bar: it == 1
                        baz: it >1
                    }
        ''')
        then:
        result.size() == 3
        and:
        result[0].val == 0
        result[0].val == Channel.STOP
        and:
        result[1].val == 1
        result[1].val == Channel.STOP
        and:
        result[2].val == 2
        result[2].val == Channel.STOP
    }

    def 'should branch and capture default' () {
        when:
        def result = dsl_eval('''   
            Channel
                .from(10,20,30)
                .branch {
                        foo: it <=10
                        bar: true
                    }
        ''')
        then:
        result.size() == 2
        and:
        result[0].val == 10
        result[0].val == Channel.STOP
        and:
        result[1].val == 20
        result[1].val == 30
        result[1].val == Channel.STOP

    }

    def 'should branch and return empty channel' () {
        when:
        def result = dsl_eval('''   
            Channel
                .from(1,2,3)
                .branch {
                        foo: it <=10
                        bar: true
                    }
        ''')
        then:
        result.size() == 2
        and:
        result[0].val == 1
        result[0].val == 2
        result[0].val == 3
        result[0].val == Channel.STOP
        and:
        result[1].val == Channel.STOP

    }


    def 'should branch and set' () {
        when:
        dsl_eval('''   
            Channel
                .from(1,2,3,40,50)
                .branch {
                    small: it < 10
                    large: it > 10
                }
                .set { result }
                
             result.small.view { "small:$it" }
             result.large.view { "large:$it" }
        ''')
        def stdout = capture.toString()
        then:
        stdout.contains('small:1')
        stdout.contains('small:2')
        stdout.contains('small:3')
        and:
        stdout.contains('large:40')
        stdout.contains('large:50')
    }

    def 'should branch and pipe set' () {
        when:
        dsl_eval('''
            Channel
                .from(1,2,3,40,50) \
                | branch { small: it < 10; large: it > 10 } \
                | set { result }
                
             result.small.view { "small:$it" }
             result.large.view { "large:$it" }
        ''')
        def stdout = capture.toString()
        then:
        stdout.contains('small:1')
        stdout.contains('small:2')
        stdout.contains('small:3')
        and:
        stdout.contains('large:40')
        stdout.contains('large:50')
    }


    def 'should branch and return custom values' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(0,1,2)
                .branch {
                        foo: it <1
                        bar: it == 1; return 10
                        baz: it >1; return 20
                    }
        ''')
        then:
        result.size() == 3
        and:
        result[0].val == 0
        result[0].val == Channel.STOP
        and:
        result[1].val == 10
        result[1].val == Channel.STOP
        and:
        result[2].val == 20
        result[2].val == Channel.STOP
    }

    def 'should handle complex nested return statement' () {
        when:
        def result = dsl_eval('''   
            Channel
                .from(-1,0,1)
                .branch {
                        foo: true
                        if( it == 0 ) { return 'zero' } 
                        else if( it<0 ) return 'less than zero'
                        else { return 'great than zero' }
                    }
        ''')
        then:
        result.val == 'less than zero'
        result.val == 'zero'
        result.val == 'great than zero'
        result.val == Channel.STOP
    }

    @Ignore // this is not supported and require explicit use of `return` 
    def 'should handle complex expression statement' () {
        when:
        def result = dsl_eval('''   
            Channel
                .from(-1,0,1)
                .branch {
                        foo: true
                        if( it == 0 ) { 'zero' } 
                        else if( it<0 ) 'less than zero'
                        else { 'great than zero' }
                    }
        ''')
        then:
        result.val == 'less than zero'
        result.val == 'zero'
        result.val == 'great than zero'
        result.val == Channel.STOP
    }



    def 'should branch and return last expression' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(0,1,2)
                .branch {
                        foo: it <1
                        bar: it == 1; it * 2 + it
                        baz: it >1; it * 2 + it 
                    }
        ''')
        then:
        result.size() == 3
        and:
        result[0].val == 0
        result[0].val == Channel.STOP
        and:
        result[1].val == 3
        result[1].val == Channel.STOP
        and:
        result[2].val == 6
        result[2].val == Channel.STOP
    }

    def 'should branch on pair argument' () {

        when:
        def result = dsl_eval('''   
            Channel
                .from(['a', 1], ['b', 2])
                .branch { key, value -> 
                        foo: key=='a'; return value
                        bar: true
                    }
        ''')
        then:
        result.size() == 2
        and:
        result[0].val == 1
        result[0].val == Channel.STOP
        and:
        result[1].val == ['b', 2]
        result[1].val == Channel.STOP
    }

    def 'should pass criteria as argument' () {
        when:
        dsl_eval('''   
            criteria = branchCriteria { 
                foo: it<5
                bar: it>=5
            }

            bra1 = Channel.from(1,2,3).branch(criteria)  
            bra2 = Channel.from(6,7,8).branch(criteria)  
            
            bra1.foo.view { "foo:$it" }
            bra2.bar.view { "bar:$it" }
        ''')
        
        def stdout = capture.toString()
        then:
        stdout.contains('foo:1')
    }

    def 'should error since param is missing' () {
        when:
        dsl_eval('''   
            Channel
                .from(1)
                .branch { -> 
                        foo: false
                        bar: true
                    }
        ''')
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains 'Branch evaluation closure should declare at least one parameter or use the implicit `it` parameter'
    }

    def 'should error due to dup label' () {
        when:
        dsl_eval('''   
            Channel.empty() .branch { foo: true; foo: true }
        ''')
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains 'Branch label already used: foo'
    }

    def 'should error due to invalid bool expr' () {
        when:
        dsl_eval('''   
            Channel.empty() .branch { foo: if(it) {}; bar: true }
        ''')
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains 'Unexpected statement'
    }

    def 'should error due to missing branch expression' () {
        when:
        dsl_eval('''   
            Channel.empty() .branch { def x=1 }
        ''')
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains 'Branch evaluation closure should contain at least one branch expression'

        when:
        dsl_eval('''   
            Channel.empty() .branch {  }
        ''')
        then:
        def e2 = thrown(ScriptCompilationException)
        e2.message.contains 'Branch evaluation closure should contain at least one branch expression'
    }


    def 'should error due to missing expression stmt' () {
        when:
        dsl_eval('''   
            Channel.empty() .branch { foo: true; if(x) {}  }
        ''')
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains 'Unexpected statement in branch condition'

    }


    def 'should branch value ch' () {

        when:
        def result = dsl_eval('''   
            Channel
                .value(10)
                .branch {
                        foo: it <5
                        otherwise: it
                    }
        ''')
        then:
        result.size() == 2
        and:
        result[0] instanceof DataflowVariable
        result[0].val == Channel.STOP
        and:
        result[1] instanceof DataflowVariable
        result[1].val == 10

    }
}
