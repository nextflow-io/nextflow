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

package nextflow.executor

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MoabExecutorTest extends Specification {

    def 'should get submit command' () {
        given:
        def exec = new MoabExecutor()
        when:
        def cmd = exec.getSubmitCommandLine(Mock(TaskRun), Paths.get('/work/dir/script.sh'))
        then:
        cmd == ['msub', '--xml', 'script.sh']
    }

    def 'should parse job id' () {
        given:
        def exec = new MoabExecutor()
        when:
        def result = exec.parseJobId('<Data><job JobID="881218"/></Data>')
        then:
        result == '881218'
        result instanceof String

        when:
        exec.parseJobId('<Data></Data>')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Missing Moab submit job ID:\n<Data></Data>\n\n'
    }

    def 'should cancel a task' () {
        given:
        def exec = new MoabExecutor()
        
        expect:
        exec.killTaskCommand('123') == ['mjobctl', '-c', '123']
        exec.killTaskCommand(['111','222','333']) == ['mjobctl', '-c', '111,222,333']
    }

    def 'should get status command' () {
        given:
        def exec = new MoabExecutor()
        when:
        def cmd = exec.queueStatusCommand('foo')
        then:
        cmd == ['showq', '--xml', '-w', "user="+System.getProperty('user.name')]
    }


    def 'should parse queue xml' () {
        given:
        def exec = new MoabExecutor()
        def STATUS = '''\
                <?xml version="1.0" encoding="UTF-8"?>
                <Data>
                   <Object>queue</Object>
                   <cluster LocalActiveNodes="479" LocalAllocProcs="16" LocalConfigNodes="747" LocalIdleNodes="229" LocalIdleProcs="3910" LocalUpNodes="708" LocalUpProcs="11968" RemoteActiveNodes="0" RemoteAllocProcs="0" RemoteConfigNodes="0" RemoteIdleNodes="0" RemoteIdleProcs="0" RemoteUpNodes="0" RemoteUpProcs="0" time="1562925910" />
                   <queue count="1" option="active">
                      <job AWDuration="1392" Account="bw18k008" Class="single" DRMJID="881196.admin2" GJID="881196" Group="hd_hd" JobID="881196" JobName="JOBNAME.deeptools_bamCompare.cell=HMG007,treat=repl1,chip=PU1,scale=SES,ratio=subtract" MasterHost="m13s0703" PAL="torque" ReqAWDuration="10800" ReqNodes="1" ReqProcs="16" RsvStartTime="1562924509" RunPriority="99" StartPriority="99" StartTime="1562924509" StatPSDed="22271.200000" StatPSUtl="1829.560100" State="Running" SubmissionTime="1562924505" SuspendDuration="0" User="hd_ow424" />
                      <job AWDuration="1392" Account="bw18k008" Class="single" DRMJID="881197.admin2" GJID="881196" Group="hd_hd" JobID="881197" JobName="JOBNAME.deeptools_bamCompare.cell=HMG007,treat=repl1,chip=PU1,scale=SES,ratio=subtract" MasterHost="m13s0703" PAL="torque" ReqAWDuration="10800" ReqNodes="1" ReqProcs="16" RsvStartTime="1562924509" RunPriority="99" StartPriority="99" StartTime="1562924509" StatPSDed="22271.200000" StatPSUtl="1829.560100" State="Running" SubmissionTime="1562924505" SuspendDuration="0" User="hd_ow424" />
                   </queue>
                   <queue count="1" option="eligible" >
                      <job AWDuration="1392" Account="bw18k008" Class="single" DRMJID="881198.admin2" GJID="881196" Group="hd_hd" JobID="881198" JobName="JOBNAME.deeptools_bamCompare.cell=HMG007,treat=repl1,chip=PU1,scale=SES,ratio=subtract" MasterHost="m13s0703" PAL="torque" ReqAWDuration="10800" ReqNodes="1" ReqProcs="16" RsvStartTime="1562924509" RunPriority="99" StartPriority="99" StartTime="1562924509" StatPSDed="22271.200000" StatPSUtl="1829.560100" State="Running" SubmissionTime="1562924505" SuspendDuration="0" User="hd_ow424" />
                   </queue>
                   <queue count="1" option="blocked" >
                      <job AWDuration="1392" Account="bw18k008" Class="single" DRMJID="881200.admin2" GJID="881200" Group="hd_hd" JobID="881200" JobName="JOBNAME.deeptools_bamCompare.cell=HMG007,treat=repl1,chip=PU1,scale=SES,ratio=subtract" MasterHost="m13s0703" PAL="torque" ReqAWDuration="10800" ReqNodes="1" ReqProcs="16" RsvStartTime="1562924509" RunPriority="99" StartPriority="99" StartTime="1562924509" StatPSDed="22271.200000" StatPSUtl="1829.560100" State="BatchHold" SubmissionTime="1562924505" SuspendDuration="0" User="hd_ow424" />
                   </queue>
                </Data>
                '''
                .stripIndent()

        when:
        def state = exec.parseQueueStatus(STATUS)
        then:
        state['881196'] == AbstractGridExecutor.QueueStatus.RUNNING
        state['881197'] == AbstractGridExecutor.QueueStatus.RUNNING
        state['881198'] == AbstractGridExecutor.QueueStatus.PENDING
        state['881200'] == AbstractGridExecutor.QueueStatus.ERROR
    }


    def "should validate headers"() {

        setup:
        def executor = [:] as MoabExecutor

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir')
        task.name = 'task name'

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                #MSUB -q alpha
                #MSUB -l walltime=00:01:00
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.time = '1m'
        task.config.memory = '1m'
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                #MSUB -q alpha
                #MSUB -l walltime=00:01:00
                #MSUB -l mem=1mb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()



        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '10m'
        task.config.memory = '5m'
        task.config.cpus = 2
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                #MSUB -q delta
                #MSUB -l nodes=1:ppn=2
                #MSUB -l walltime=00:10:00
                #MSUB -l mem=5mb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '1d'
        task.config.memory = '1g'
        task.config.cpus = 8
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                #MSUB -q delta
                #MSUB -l nodes=1:ppn=8
                #MSUB -l walltime=24:00:00
                #MSUB -l mem=1gb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '2d 6h 10m'
        task.config.memory = '2g'
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o /work/dir/.command.log
                #MSUB -j oe
                #MSUB -q delta
                #MSUB -l walltime=54:10:00
                #MSUB -l mem=2gb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

    }

    def WorkDirWithBlanks() {

        setup:
        def executor = Spy(MoabExecutor)

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir 1')
        task.name = 'task name'

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #MSUB -N nf-task_name
                #MSUB -o "/work/dir\\ 1/.command.log"
                #MSUB -j oe
                NXF_CHDIR=/work/dir\\ 1
                '''
                .stripIndent().leftTrim()

    }
}
