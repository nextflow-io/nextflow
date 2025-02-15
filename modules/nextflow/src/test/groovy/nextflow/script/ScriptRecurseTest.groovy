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
 *
 */

package nextflow.script

import nextflow.Channel
import nextflow.NextflowMeta
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class ScriptRecurseTest extends Dsl2Spec {

    def setupSpec() { NextflowMeta.instance.preview.recursion=true }
    def cleanupSpec() { NextflowMeta.instance.preview.recursion=false }

    def 'should recourse a process execution' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: val x
          output: val y
          exec: y = x +1
        }
       
        workflow {
            main: foo.recurse(1).times(3)
            emit: foo.out
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 2
        result.val == 3
        result.val == 4
        result.val == Channel.STOP
    }

    def 'should recourse a process until a condition is verified' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: val x
          output: val y
          exec: y = x+1
        }
       
        workflow {
            main: foo.recurse(1).until { it >= 4 }
            emit: foo.out
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 2
        result.val == 3
        result.val == 4
        result.val == Channel.STOP
    }


    def 'should recourse a workflow execution' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: val x
          output: val y
          exec: y = x+1
        }
  
        process bar {
          input: val x
          output: val y
          exec: y = x*x
        }
       
        workflow group {
            take: x 
            main: 
              foo(x)
              bar(foo.out)
            emit:
              bar.out
        }
        
        workflow {
            main: group.recurse(1).times(3)
            emit: group.out  
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 4
        result.val == 25
        result.val == 676
        result.val == Channel.STOP
    }

    def 'should recurse with scan' () {
        given:
        def SCRIPT = '''
         process foo {
           input:
             val x
           output:
             val z
           exec:
              z = x.sum()+1
         }
         
         workflow {
           main:
             data = channel.of(10,20,30) 
             foo.scan(data)
           emit:
             foo.out
         }
        '''
        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 11 // 10 +1
        result.val == 32 // 20 + 11 +1
        result.val == 74 // 30 + 11 + 32 +1

    }

}
