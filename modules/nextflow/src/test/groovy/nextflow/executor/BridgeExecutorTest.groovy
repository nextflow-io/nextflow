/*
 * Copyright 2013-2024, Seqera Labs
 * Copyright 2022, CEA-CNRGH
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
import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Eric Bonnet <eric.d.bonnet@gmail.com>
 */
class BridgeExecutorTest extends Specification {

    def testParseJob() {

        given:
        def exec = [:] as BridgeExecutor 

        expect:
        //submission pattern example: Submitted Batch Session 1277017
        exec.parseJobId('Submitted Batch Session 10') == '10'
        exec.parseJobId('Submitted Batch Session 20') == '20'

        when:
        exec.parseJobId('Something else 10')
        then:
        thrown(IllegalStateException)

    }

    def testKill() {

        given:
        def exec = [:] as BridgeExecutor 
        expect:
        exec.killTaskCommand(123) == ['ccc_mdel','123']

    }

    def testGetCommandLine() {

        when:
        def exec = [:] as BridgeExecutor 
        then:
        exec.getSubmitCommandLine(Mock(TaskRun), Paths.get('/some/path/job.sh')) == ['ccc_msub', 'job.sh']
    }

    def testGetHeaders() {

        setup:
        // Bridge executor
        def executor = [:] as BridgeExecutor

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/path')
        task.name = 'virtual_digest'

        when:
        task.index = 21
        task.config = new TaskConfig()
        task.config.queue = 'normal'
        task.config.time = '10h'
        task.config.memory = '8 GB' 
        task.config.cpus = 2
        then:
        executor.getHeaders(task) == '''
                #MSUB -r nf-virtual_digest
                #MSUB -o /work/path/.command.log
                #MSUB -c 2
                #MSUB -T 36000
                #MSUB -M 8192
                #MSUB -q normal
                '''
                .stripIndent().leftTrim()
    }

    def testQstatCommand() {

        setup:
        def executor = [:] as BridgeExecutor
        def text =
                """
                5 pending
                6 pending
                13 running
                15 failed
                4 running
                22 suspended
                7 unknown
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 7
        result['4'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5'] == AbstractGridExecutor.QueueStatus.PENDING
        result['6'] == AbstractGridExecutor.QueueStatus.PENDING
        result['13'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['15'] == AbstractGridExecutor.QueueStatus.ERROR
        result['7'] == AbstractGridExecutor.QueueStatus.ERROR
        result['22'] == AbstractGridExecutor.QueueStatus.HOLD

    }

}
