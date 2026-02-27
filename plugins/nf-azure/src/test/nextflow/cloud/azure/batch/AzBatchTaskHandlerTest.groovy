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

package nextflow.cloud.azure.batch

import com.azure.compute.batch.models.BatchTask
import com.azure.compute.batch.models.BatchTaskExecutionInfo
import com.azure.compute.batch.models.BatchTaskFailureInfo
import com.azure.compute.batch.models.BatchTaskState
import com.azure.compute.batch.models.ErrorCategory
import com.sun.jna.platform.unix.X11
import nextflow.processor.TaskStatus

import java.nio.file.Path

import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.Executor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchTaskHandlerTest extends Specification {

    def createTaskRun() {
        Mock(TaskRun) {
            name >> 'foo'
            workDir >> Path.of('/work/dir')
            container >> 'ubuntu'
        }
    }

    def 'should validate config' () {
        given:
        def exec = Mock(AzBatchExecutor)

        when:
        def task = Mock(TaskRun) {
            name >> 'foo'
            workDir >> Path.of('/work/dir')
        }
        and:
        new AzBatchTaskHandler(task, exec)
        then:
        def e = thrown(ProcessUnrecoverableException)
        e.message.startsWith('No container image specified for process foo')

        when:
        task = createTaskRun()
        and:
        new AzBatchTaskHandler(task, exec)
        then:
        noExceptionThrown()
    }

    def 'should submit task'() {
        given:
        def azure = Mock(AzBatchService)
        def executor = Mock(AzBatchExecutor)
        def processor = Mock(TaskProcessor) {
            getExecutor() >> executor
        }
        def task = createTaskRun()
        task.getProcessor() >> processor
        task.getConfig() >> Mock(TaskConfig)
        and:
        def handler = Spy(new AzBatchTaskHandler(task, executor)) {
            getBatchService() >> azure
        }

        when:
        handler.submit()

        then:
        1 * handler.createBashWrapper() >> Mock(BashWrapperBuilder)
        1 * handler.getBatchService() >> Mock(AzBatchService)
    }

    def 'should create the trace record' () {
        given:
        def exec = Mock(AzBatchExecutor) { getName() >> 'azurebatch' }
        def processor = Mock(TaskProcessor)
        processor.getExecutor() >> exec
        processor.getName() >> 'foo'
        processor.getConfig() >> new ProcessConfig(Mock(BaseScript))
        def task = createTaskRun()
        task.getProcessor() >> processor
        task.getConfig() >> GroovyMock(TaskConfig)
        def handler = Spy(new AzBatchTaskHandler(task, exec))
        handler.@taskKey = new AzTaskKey('job-123', 'nf-456')

        when:
        def trace = handler.getTraceRecord()
        then:
        1 * handler.isCompleted() >> false
        1 * handler.getMachineInfo() >> new CloudMachineInfo('Standard1', 'west-eu', PriceModel.standard)

        and:
        trace.native_id == 'job-123/nf-456'
        trace.executorName == 'azurebatch'
        trace.machineInfo.type == 'Standard1'
        trace.machineInfo.zone == 'west-eu'
        trace.machineInfo.priceModel == PriceModel.standard
    }

    def 'should check if completed with exit code from scheduler'() {
        given:
        def task = Spy(new TaskRun()){
            getContainer() >> 'ubuntu'
        }
        task.name = 'foo'
        task.workDir = Path.of('/tmp/wdir')
        def taskKey = new AzTaskKey('pool-123', 'job-456')
        def azTask = new BatchTask()
        def execInfo = new BatchTaskExecutionInfo(0,0)
        execInfo.exitCode = 0
        azTask.executionInfo = execInfo
        azTask.state = BatchTaskState.COMPLETED

        def batchService = Mock(AzBatchService){
            getTask(taskKey) >> azTask
        }
        def executor = Mock(AzBatchExecutor){
            getBatchService() >> batchService
        }
        def handler = Spy(new AzBatchTaskHandler(task, executor)){
            deleteTask(_,_) >> null
        }
        handler.status = TaskStatus.RUNNING
        handler.taskKey = taskKey

        when:
        def result = handler.checkIfCompleted()
        then:
        0 * handler.readExitFile()  // Should NOT read exit file when scheduler provides exit code
        and:
        result == true
        handler.task.exitStatus == 0
        handler.status == TaskStatus.COMPLETED

    }

    def 'should check if completed with non-zero exit code from scheduler'() {
        given:
        def task = Spy(new TaskRun()){
            getContainer() >> 'ubuntu'
        }
        task.name = 'foo'
        task.workDir = Path.of('/tmp/wdir')
        def taskKey = new AzTaskKey('pool-123', 'job-456')
        def azTask = new BatchTask()
        def execInfo = new BatchTaskExecutionInfo(0,0)
        execInfo.exitCode = 137
        azTask.executionInfo = execInfo
        azTask.state = BatchTaskState.COMPLETED

        def batchService = Mock(AzBatchService){
            getTask(taskKey) >> azTask
        }
        def executor = Mock(AzBatchExecutor){
            getBatchService() >> batchService
        }
        def handler = Spy(new AzBatchTaskHandler(task, executor)){
            deleteTask(_,_) >> null
        }
        handler.status = TaskStatus.RUNNING
        handler.taskKey = taskKey

        when:
        def result = handler.checkIfCompleted()
        then:
        0 * handler.readExitFile()  // Should NOT read exit file when scheduler provides exit code
        and:
        result == true
        handler.task.exitStatus == 137
        handler.status == TaskStatus.COMPLETED


    }

    def 'should check if completed and fallback to exit file when scheduler exit code is null'() {
        given:
        def task = Spy(new TaskRun()){
            getContainer() >> 'ubuntu'
        }
        task.name = 'foo'
        task.workDir = Path.of('/tmp/wdir')
        def taskKey = new AzTaskKey('pool-123', 'job-456')
        def azTask = new BatchTask()
        def execInfo = new BatchTaskExecutionInfo(0,0)
        azTask.executionInfo = execInfo
        azTask.state = BatchTaskState.COMPLETED

        def batchService = Mock(AzBatchService){
            getTask(taskKey) >> azTask
        }
        def executor = Mock(AzBatchExecutor){
            getBatchService() >> batchService
        }
        def handler = Spy(new AzBatchTaskHandler(task, executor)){
            deleteTask(_,_) >> null
        }
        handler.status = TaskStatus.RUNNING
        handler.taskKey = taskKey

        when:
        def result = handler.checkIfCompleted()

        then:
        1 * handler.readExitFile() >> 0     // Should read exit file as fallback
        and:
        result == true
        handler.task.exitStatus == 0
        handler.status == TaskStatus.COMPLETED
    }

    def 'should check if completed and no scheduler exit code neither .exitcode file'() {
        given:
        def task = Spy(new TaskRun()){
            getContainer() >> 'ubuntu'
        }
        task.name = 'foo'
        task.workDir = Path.of('/tmp/wdir')
        def taskKey = new AzTaskKey('pool-123', 'job-456')
        def azTask = new BatchTask()
        def execInfo = new BatchTaskExecutionInfo(0,0)
        def failureInfo = new BatchTaskFailureInfo(ErrorCategory.USER_ERROR)
        failureInfo.message = 'Unknown error'
        execInfo.failureInfo = failureInfo
        azTask.executionInfo = execInfo
        azTask.state = BatchTaskState.COMPLETED

        def batchService = Mock(AzBatchService){
            getTask(taskKey) >> azTask
        }
        def executor = Mock(AzBatchExecutor){
            getBatchService() >> batchService
        }
        def handler = Spy(new AzBatchTaskHandler(task, executor)){
            deleteTask(_,_) >> null
        }
        handler.status = TaskStatus.RUNNING
        handler.taskKey = taskKey

        when:
        def result = handler.checkIfCompleted()

        then:
        1 * handler.readExitFile() >> Integer.MAX_VALUE     // Should read exit file as fallback
        and:
        result == true
        handler.task.exitStatus == Integer.MAX_VALUE
        handler.status == TaskStatus.COMPLETED
        handler.task.error.message == 'Unknown error'
    }
}
