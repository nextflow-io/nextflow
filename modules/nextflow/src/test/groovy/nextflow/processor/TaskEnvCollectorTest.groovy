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

package nextflow.processor

import nextflow.exception.ProcessEvalException
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskEnvCollectorTest extends Specification {

    def 'should parse env map' () {
        given:
        def workDir = TestHelper.createInMemTempDir()
        def envFile = workDir.resolve(TaskRun.CMD_ENV)
        envFile.text =  '''
            ALPHA=one
            /ALPHA/
            DELTA=x=y
            /DELTA/
            OMEGA=
            /OMEGA/
            LONG=one
            two
            three
            /LONG/=exit:0
            '''.stripIndent()

        when:
        def result = new TaskEnvCollector(workDir, Map.of()).collect()
        then:
        result == [ALPHA:'one', DELTA: "x=y", OMEGA: '', LONG: 'one\ntwo\nthree']
    }

    def 'should parse env map with command error' () {
        given:
        def workDir = TestHelper.createInMemTempDir()
        def envFile = workDir.resolve(TaskRun.CMD_ENV)
        envFile.text =  '''
            ALPHA=one
            /ALPHA/
            cmd_out_1=Hola
            /cmd_out_1/=exit:0
            cmd_out_2=This is an error message
            for unknown reason
            /cmd_out_2/=exit:100
            '''.stripIndent()

        when:
        new TaskEnvCollector(workDir, [cmd_out_1: 'foo --this', cmd_out_2: 'bar --that']).collect()
        then:
        def e = thrown(ProcessEvalException)
        e.message == 'Unable to evaluate output'
        e.command == 'bar --that'
        e.output == 'This is an error message\nfor unknown reason'
        e.status == 100
    }

}
