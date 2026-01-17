/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.executor

import io.seqera.sched.api.schema.v1a1.MachineInfo
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.TaskAttempt
import io.seqera.sched.api.schema.v1a1.TaskState as SchedTaskState
import io.seqera.sched.api.schema.v1a1.TaskStatus as SchedTaskStatus
import io.seqera.sched.client.SchedClient
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import spock.lang.Specification
import spock.lang.Unroll

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

    @Unroll
    def 'should convert price model #schedModel to #expected'() {
        given:
        def handler = createHandler()

        expect:
        handler.convertPriceModel(schedModel) == expected

        where:
        schedModel              | expected
        SchedPriceModel.SPOT    | PriceModel.spot
        SchedPriceModel.STANDARD| PriceModel.standard
        null                    | null
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
}
