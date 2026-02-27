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

package io.seqera.executor

import com.google.common.hash.HashCode
import io.seqera.config.ExecutorOpts
import io.seqera.sched.api.schema.v1a1.DescribeTaskResponse
import io.seqera.sched.api.schema.v1a1.GetTaskLogsResponse
import io.seqera.sched.api.schema.v1a1.MachineInfo
import io.seqera.sched.api.schema.v1a1.NextflowTask
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ResourceLimit
import io.seqera.sched.api.schema.v1a1.ResourceRequirement
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.api.schema.v1a1.TaskAttempt
import io.seqera.sched.api.schema.v1a1.TaskState as SchedTaskState
import io.seqera.sched.api.schema.v1a1.TaskStatus as SchedTaskStatus
import io.seqera.sched.client.SchedClient
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.exception.ProcessException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import spock.lang.Specification

import java.nio.file.Paths

/**
 * Tests for SeqeraTaskHandler metadata fetching functionality
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SeqeraTaskHandlerTest extends Specification {

    def 'should return null for getMachineInfo when cachedTaskState is null'() {
        given:
        def handler = createHandler()

        expect:
        handler.getMachineInfo() == null
    }

    def 'should return null for getMachineInfo when attempts list is empty'() {
        given:
        def handler = createHandler()
        handler.cachedTaskState = new SchedTaskState().attempts([])

        expect:
        handler.getMachineInfo() == null
    }

    def 'should return null for getMachineInfo when last attempt has no machine info'() {
        given:
        def handler = createHandler()
        def attempt = new TaskAttempt()
            .index(1)
            .nativeId('arn:aws:ecs:us-east-1:123:task/abc')
            .status(SchedTaskStatus.SUCCEEDED)
        handler.cachedTaskState = new SchedTaskState().attempts([attempt])

        expect:
        handler.getMachineInfo() == null
    }

    def 'should extract machine info from last task attempt'() {
        given:
        def handler = createHandler()
        def machineInfo = new MachineInfo()
            .type('m5.large')
            .zone('us-east-1a')
            .priceModel(SchedPriceModel.SPOT)
        def attempt = new TaskAttempt()
            .index(1)
            .nativeId('arn:aws:ecs:us-east-1:123:task/abc')
            .status(SchedTaskStatus.SUCCEEDED)
            .machineInfo(machineInfo)
        handler.cachedTaskState = new SchedTaskState().attempts([attempt])

        when:
        def result = handler.getMachineInfo()

        then:
        result != null
        result.type == 'm5.large'
        result.zone == 'us-east-1a'
        result.priceModel == PriceModel.spot
    }

    def 'should cache machine info after first retrieval'() {
        given:
        def handler = createHandler()
        def machineInfo = new MachineInfo()
            .type('c5.xlarge')
            .zone('eu-west-1b')
            .priceModel(SchedPriceModel.STANDARD)
        def attempt = new TaskAttempt()
            .index(1)
            .nativeId('arn:aws:ecs:eu-west-1:123:task/xyz')
            .status(SchedTaskStatus.SUCCEEDED)
            .machineInfo(machineInfo)
        handler.cachedTaskState = new SchedTaskState().attempts([attempt])

        when:
        def first = handler.getMachineInfo()
        // Clear the cached task state
        handler.cachedTaskState = null
        def second = handler.getMachineInfo()

        then:
        first.is(second)
    }

    def 'should use last attempt when multiple attempts exist'() {
        given:
        def handler = createHandler()
        def info1 = new MachineInfo().type('m5.large').zone('us-east-1a').priceModel(SchedPriceModel.SPOT)
        def info2 = new MachineInfo().type('c5.xlarge').zone('us-east-1b').priceModel(SchedPriceModel.STANDARD)
        def attempt1 = new TaskAttempt().index(1).nativeId('task-1').status(SchedTaskStatus.FAILED).machineInfo(info1)
        def attempt2 = new TaskAttempt().index(2).nativeId('task-2').status(SchedTaskStatus.SUCCEEDED).machineInfo(info2)
        handler.cachedTaskState = new SchedTaskState().attempts([attempt1, attempt2])

        when:
        def result = handler.getMachineInfo()

        then:
        result.type == 'c5.xlarge'
        result.zone == 'us-east-1b'
        result.priceModel == PriceModel.standard
    }

    def 'should return null for getNumSpotInterruptions when task not completed'() {
        given:
        def handler = createHandler()
        handler.taskId = 'task-123'
        handler.status = TaskStatus.RUNNING
        handler.cachedTaskState = new SchedTaskState().numSpotInterruptions(2)

        expect:
        handler.getNumSpotInterruptions() == null
    }

    def 'should return null for getNumSpotInterruptions when cachedTaskState is null'() {
        given:
        def handler = createHandler()
        handler.taskId = 'task-123'
        handler.status = TaskStatus.COMPLETED

        expect:
        handler.getNumSpotInterruptions() == null
    }

    def 'should return spot interruptions count from task state'() {
        given:
        def handler = createHandler()
        handler.taskId = 'task-123'
        handler.status = TaskStatus.COMPLETED
        handler.cachedTaskState = new SchedTaskState().numSpotInterruptions(3)

        expect:
        handler.getNumSpotInterruptions() == 3
    }

    def 'should return zero spot interruptions when none occurred'() {
        given:
        def handler = createHandler()
        handler.taskId = 'task-123'
        handler.status = TaskStatus.COMPLETED
        handler.cachedTaskState = new SchedTaskState().numSpotInterruptions(0)

        expect:
        handler.getNumSpotInterruptions() == 0
    }

    def 'should return null for getLogStreamId when cachedTaskState is null'() {
        given:
        def handler = createHandler()

        expect:
        handler.getLogStreamId() == null
    }

    def 'should return log stream id from cached task state'() {
        given:
        def handler = createHandler()
        handler.cachedTaskState = new SchedTaskState().logStreamId('log-stream-abc123')

        expect:
        handler.getLogStreamId() == 'log-stream-abc123'
    }

    def 'should return null for getNativeId when cachedTaskState is null'() {
        given:
        def handler = createHandler()

        expect:
        handler.getNativeId() == null
    }

    def 'should return task id from cached task state'() {
        given:
        def handler = createHandler()
        handler.cachedTaskState = new SchedTaskState().id('tsk-abc123')

        expect:
        handler.getNativeId() == 'tsk-abc123'
    }

    def 'should populate trace record with metadata'() {
        given:
        def handler = createHandlerForTraceTest()
        handler.taskId = 'task-123'
        handler.status = TaskStatus.COMPLETED
        def machineInfo = new MachineInfo()
            .type('m5.large')
            .zone('us-east-1a')
            .priceModel(SchedPriceModel.SPOT)
        def attempt = new TaskAttempt()
            .index(1)
            .nativeId('arn:aws:ecs:us-east-1:123:task/abc')
            .status(SchedTaskStatus.SUCCEEDED)
            .machineInfo(machineInfo)
        handler.cachedTaskState = new SchedTaskState()
            .id('tsk-xyz789')
            .attempts([attempt])
            .numSpotInterruptions(2)
            .logStreamId('log-stream-xyz')

        when:
        def trace = handler.getTraceRecord()

        then:
        trace.get('native_id') == 'tsk-xyz789'
        trace.getMachineInfo().type == 'm5.large'
        trace.getNumSpotInterruptions() == 2
        trace.getLogStreamId() == 'log-stream-xyz'
        trace.getExecutorName() == 'seqera/aws'
    }

    def 'should detect completion when batch submission fails'() {
        given:
        def error = new RuntimeException('Batch submission failed')
        def handler = createHandlerWithError(error)
        // Simulate batch submission failure - task has error and status is COMPLETED but never went through RUNNING
        handler.status = TaskStatus.COMPLETED

        expect:
        // checkIfCompleted should return true to allow proper error propagation
        handler.checkIfCompleted() == true
    }

    def 'should not detect completion for normal pending task'() {
        given:
        def handler = createHandler()
        // Normal task that has been submitted but not yet running
        handler.status = TaskStatus.SUBMITTED

        expect:
        // checkIfCompleted should return false - task is still pending
        handler.checkIfCompleted() == false
    }

    def 'should set task error with error message when failed with no exit code'() {
        given:
        Throwable capturedError = null
        Integer capturedExitStatus = null
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> Mock(TaskConfig)
            lazyName() >> 'test_task'
            setExitStatus(_) >> { args -> capturedExitStatus = args[0] }
            getExitStatus() >> { capturedExitStatus }
            setError(_) >> { args -> capturedError = args[0] }
            setStdout(_) >> {}
            setStderr(_) >> {}
        }
        def taskState = new SchedTaskState()
            .status(SchedTaskStatus.FAILED)
            .errorMessage('Container terminated with OOMKilled')
        def describeResponse = new DescribeTaskResponse().taskState(taskState)
        def client = Mock(SchedClient) {
            describeTask(_) >> describeResponse
            getTaskLogs(_) >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> client
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            readExitFile() >> Integer.MAX_VALUE
        }
        handler.setBatchTaskId('task-123')
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then:
        completed
        capturedExitStatus == Integer.MAX_VALUE
        capturedError instanceof ProcessException
        capturedError.message == 'Container terminated with OOMKilled'
    }

    def 'should set task error with fallback message when failed with no exit code and no error message'() {
        given:
        Throwable capturedError = null
        Integer capturedExitStatus = null
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> Mock(TaskConfig)
            lazyName() >> 'test_task'
            setExitStatus(_) >> { args -> capturedExitStatus = args[0] }
            getExitStatus() >> { capturedExitStatus }
            setError(_) >> { args -> capturedError = args[0] }
            setStdout(_) >> {}
            setStderr(_) >> {}
        }
        def taskState = new SchedTaskState()
            .status(SchedTaskStatus.FAILED)
        def describeResponse = new DescribeTaskResponse().taskState(taskState)
        def client = Mock(SchedClient) {
            describeTask(_) >> describeResponse
            getTaskLogs(_) >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> client
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            readExitFile() >> Integer.MAX_VALUE
        }
        handler.setBatchTaskId('task-456')
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then:
        completed
        capturedExitStatus == Integer.MAX_VALUE
        capturedError instanceof ProcessException
        capturedError.message == 'Task failed for unknown reason'
    }

    def 'should not set task error when failed with valid exit code'() {
        given:
        Throwable capturedError = null
        Integer capturedExitStatus = null
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> Mock(TaskConfig)
            lazyName() >> 'test_task'
            setExitStatus(_) >> { args -> capturedExitStatus = args[0] }
            getExitStatus() >> { capturedExitStatus }
            setError(_) >> { args -> capturedError = args[0] }
            setStdout(_) >> {}
            setStderr(_) >> {}
        }
        def taskState = new SchedTaskState()
            .status(SchedTaskStatus.FAILED)
            .errorMessage('Some error')
        def describeResponse = new DescribeTaskResponse().taskState(taskState)
        def client = Mock(SchedClient) {
            describeTask(_) >> describeResponse
            getTaskLogs(_) >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> client
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            readExitFile() >> 1
        }
        handler.setBatchTaskId('task-789')
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then:
        completed
        capturedExitStatus == 1
        capturedError == null
    }

    def 'should set index and hash on submitted task'() {
        given:
        Task capturedTask = null
        def batchSubmitter = Mock(SeqeraBatchSubmitter) {
            submit(_, _) >> { handler, task -> capturedTask = task }
        }
        def seqeraConfig = Mock(ExecutorOpts) {
            getMachineRequirement() >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getBatchSubmitter() >> batchSubmitter
            getSeqeraConfig() >> seqeraConfig
        }
        def taskConfig = Mock(TaskConfig) {
            getCpus() >> 1
            getMemory() >> null
            getAccelerator() >> null
            getDisk() >> null
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> taskConfig
            getContainer() >> 'ubuntu:latest'
            getContainerPlatform() >> null
            lazyName() >> 'test_task'
            getId() >> TaskId.of(42)
            getHash() >> HashCode.fromString('abcd1234')
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionEnabled() >> true
            fusionSubmitCli() >> ['bash', '-c', 'echo hello']
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [:]
            }
        }

        when:
        handler.submit()

        then:
        1 * executor.ensureRunCreated()
        capturedTask != null
        capturedTask.getNextflow() != null
        capturedTask.getNextflow().getTaskId() == 42
        capturedTask.getNextflow().getHash() == 'abcd1234'
        capturedTask.getNextflow().getWorkDir() == '/work/ab/cd1234'
    }

    def 'should handle null task id and hash'() {
        given:
        Task capturedTask = null
        def batchSubmitter = Mock(SeqeraBatchSubmitter) {
            submit(_, _) >> { handler, task -> capturedTask = task }
        }
        def seqeraConfig = Mock(ExecutorOpts) {
            getMachineRequirement() >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getBatchSubmitter() >> batchSubmitter
            getSeqeraConfig() >> seqeraConfig
        }
        def taskConfig = Mock(TaskConfig) {
            getCpus() >> 1
            getMemory() >> null
            getAccelerator() >> null
            getDisk() >> null
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> taskConfig
            getContainer() >> 'ubuntu:latest'
            getContainerPlatform() >> null
            lazyName() >> 'test_task'
            getId() >> null
            getHash() >> null
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionEnabled() >> true
            fusionSubmitCli() >> ['bash', '-c', 'echo hello']
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [:]
            }
        }

        when:
        handler.submit()

        then:
        capturedTask != null
        capturedTask.getNextflow() != null
        capturedTask.getNextflow().getTaskId() == null
        capturedTask.getNextflow().getHash() == null
        capturedTask.getNextflow().getWorkDir() == '/work/ab/cd1234'
    }

    def 'should return only fusion env when no config environment'() {
        given:
        def seqeraConfig = Mock(ExecutorOpts) {
            getTaskEnvironment() >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getSeqeraConfig() >> seqeraConfig
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig)
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [FUSION_KEY: 'fusion_val']
            }
        }

        when:
        def result = handler.getTaskEnvironment()

        then:
        result == [FUSION_KEY: 'fusion_val']
    }

    def 'should merge config environment with fusion env'() {
        given:
        def seqeraConfig = Mock(ExecutorOpts) {
            getTaskEnvironment() >> [MY_VAR: 'my_val', OTHER: 'other_val']
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getSeqeraConfig() >> seqeraConfig
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig)
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [FUSION_KEY: 'fusion_val']
            }
        }

        when:
        def result = handler.getTaskEnvironment()

        then:
        result == [MY_VAR: 'my_val', OTHER: 'other_val', FUSION_KEY: 'fusion_val']
    }

    def 'should give fusion env precedence over config environment'() {
        given:
        def seqeraConfig = Mock(ExecutorOpts) {
            getTaskEnvironment() >> [SHARED_KEY: 'config_val', MY_VAR: 'my_val']
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getSeqeraConfig() >> seqeraConfig
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig)
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [SHARED_KEY: 'fusion_val']
            }
        }

        when:
        def result = handler.getTaskEnvironment()

        then:
        result == [MY_VAR: 'my_val', SHARED_KEY: 'fusion_val']
    }

    def 'should return granted time from resource requirement'() {
        given:
        def handler = createHandler()
        handler.cachedTaskState = new SchedTaskState()
            .resourceRequirement(new ResourceRequirement().time('2h'))

        expect:
        handler.getGrantedTime() == Duration.of('2h').toMillis()
    }

    def 'should fallback to config time when cachedTaskState is null'() {
        given:
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig) { getTime() >> Duration.of('6h') }
        }
        def executor = Mock(SeqeraExecutor) { getClient() >> Mock(SchedClient) }
        def handler = new SeqeraTaskHandler(taskRun, executor)

        expect:
        handler.getGrantedTime() == Duration.of('6h').toMillis()
    }

    def 'should build resource limit from task config'() {
        given:
        def taskConfig = Mock(TaskConfig) {
            getResourceLimit('memory') >> MemoryUnit.of('2 GB')
            getResourceLimit('cpus') >> 4
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> taskConfig
        }
        def executor = Mock(SeqeraExecutor) { getClient() >> Mock(SchedClient) }
        def handler = new SeqeraTaskHandler(taskRun, executor)

        when:
        def result = handler.toResourceLimit()

        then:
        result != null
        result.memoryMiB == 2048
        result.cpuShares == 4096
    }

    def 'should build resource limit with only memory'() {
        given:
        def taskConfig = Mock(TaskConfig) {
            getResourceLimit('memory') >> MemoryUnit.of('800 MB')
            getResourceLimit('cpus') >> null
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> taskConfig
        }
        def executor = Mock(SeqeraExecutor) { getClient() >> Mock(SchedClient) }
        def handler = new SeqeraTaskHandler(taskRun, executor)

        when:
        def result = handler.toResourceLimit()

        then:
        result != null
        result.memoryMiB == 800
        result.cpuShares == null
    }

    def 'should return null resource limit when no limits defined'() {
        given:
        def taskConfig = Mock(TaskConfig) {
            getResourceLimit('memory') >> null
            getResourceLimit('cpus') >> null
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> taskConfig
        }
        def executor = Mock(SeqeraExecutor) { getClient() >> Mock(SchedClient) }
        def handler = new SeqeraTaskHandler(taskRun, executor)

        when:
        def result = handler.toResourceLimit()

        then:
        result == null
    }

    def 'should include resource limit in submitted task'() {
        given:
        Task capturedTask = null
        def batchSubmitter = Mock(SeqeraBatchSubmitter) {
            submit(_, _) >> { handler, task -> capturedTask = task }
        }
        def seqeraConfig = Mock(ExecutorOpts) {
            getMachineRequirement() >> null
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getBatchSubmitter() >> batchSubmitter
            getSeqeraConfig() >> seqeraConfig
        }
        def taskConfig = Mock(TaskConfig) {
            getCpus() >> 1
            getMemory() >> null
            getAccelerator() >> null
            getDisk() >> null
            getResourceLimit('memory') >> MemoryUnit.of('2 GB')
            getResourceLimit('cpus') >> null
        }
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getConfig() >> taskConfig
            getContainer() >> 'ubuntu:latest'
            getContainerPlatform() >> null
            lazyName() >> 'test_task'
            getId() >> TaskId.of(1)
            getHash() >> HashCode.fromString('abcd1234')
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionEnabled() >> true
            fusionSubmitCli() >> ['bash', '-c', 'echo hello']
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [:]
            }
        }

        when:
        handler.submit()

        then:
        capturedTask != null
        capturedTask.getResourceLimit() != null
        capturedTask.getResourceLimit().memoryMiB == 2048
        capturedTask.getResourceLimit().cpuShares == null
    }

    /**
     * Creates a test handler with minimal mocked dependencies
     */
    private SeqeraTaskHandler createHandler() {
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig)
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
        }
        return new SeqeraTaskHandler(taskRun, executor)
    }

    /**
     * Creates a test handler with an error set on the task
     */
    private SeqeraTaskHandler createHandlerWithError(Throwable error) {
        def taskRun = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getConfig() >> Mock(TaskConfig)
            getError() >> error
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
        }
        return new SeqeraTaskHandler(taskRun, executor)
    }

    /**
     * Creates a test handler with mocks sufficient for getTraceRecord()
     */
    private SeqeraTaskHandler createHandlerForTraceTest() {
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getName() >> 'seqera'
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'test_process'
            getExecutor() >> executor
        }
        def taskConfig = new TaskConfig(attempt: 1, cpus: 1)
        def taskRun = Mock(TaskRun) {
            getId() >> TaskId.of(1)
            getHashLog() >> 'ab/cd1234'
            getName() >> 'test_task'
            getExitStatus() >> 0
            getProcessor() >> processor
            getConfig() >> taskConfig
            getContainer() >> 'ubuntu:latest'
            getTraceScript() >> 'echo hello'
            getScratch() >> null
            getWorkDirStr() >> '/work/ab/cd1234'
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getEnvironmentStr() >> ''
            containerMeta() >> null
        }
        return new SeqeraTaskHandler(taskRun, executor)
    }
}
