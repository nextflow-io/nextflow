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

package nextflow.script

import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.NF
import nextflow.NextflowMeta
import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Timeout
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Ignore // TODO fix these tests 
@Timeout(5)
class ScriptDsl2Test extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2(); NF.init() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }


    def 'should allow pipe process and operator' () {
        given:
        def SCRIPT = '''
        process foo {
          output: val result
          exec: result = "hello"
        }     
 
        process bar {
          output: val result
          exec: result = "world"
        } 
        
        workflow {
           emit: (foo & bar) | concat      
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'hello'
        result[0].val == 'world'
        result[0].val == Channel.STOP
    }
    
    def 'should allow process and operator composition' () {
        given:
        def SCRIPT = '''
        process foo {
          output: val result
          exec: result = "hello"
        }     
 
        process bar {
          output: val result
          exec: result = "world"
        } 
        
        workflow {
           main: foo(); bar()
           emit: foo.out.concat(bar.out)      
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'hello'
        result[0].val == 'world'
        result[0].val == Channel.STOP
    }

}
