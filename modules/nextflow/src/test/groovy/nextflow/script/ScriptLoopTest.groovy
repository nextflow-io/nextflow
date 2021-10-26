/*
 * Copyright 2020-2021, Seqera Labs
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
import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptLoopTest extends Dsl2Spec {

    def 'should define a process with output alias' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: val x
          output: val y
          exec: y = x+1
        }
       
        workflow {
            main: foo.loop(1).times(3)
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

}
