/*
 * Copyright 2025, Seqera
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

import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzBatchOpts
import nextflow.executor.Executor
import nextflow.processor.TaskProcessor
import spock.lang.Specification

/**
 * Tests for {@link AzBatchProcessObserver}
 *
 * @author Adam Talbot <adam.talbot@seqera.io>
 */
class AzBatchProcessObserverTest extends Specification {

    def 'should set job termination for matching processor'() {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(AzBatchExecutor)
        def batchService = Mock(AzBatchService)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts)
        and:
        processor.executor >> executor
        processor.name >> 'test_process'
        executor.batchService >> batchService
        batchService.config >> config
        config.batch() >> batchOpts
        batchOpts.terminateJobsOnCompletion >> true
        batchService.getJobIdsForProcessor(processor) >> ['job-1', 'job-2']

        when:
        def observer = new AzBatchProcessObserver(Mock(Session))
        observer.onProcessTerminate(processor)

        then:
        1 * batchService.setJobTermination('job-1')
        1 * batchService.setJobTermination('job-2')
    }

    def 'should skip non-AzBatch processors'() {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(Executor)
        and:
        processor.executor >> executor

        when:
        def observer = new AzBatchProcessObserver(Mock(Session))
        observer.onProcessTerminate(processor)

        then:
        noExceptionThrown()
    }

    def 'should skip when terminateJobsOnCompletion is disabled'() {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(AzBatchExecutor)
        def batchService = Mock(AzBatchService)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts)
        and:
        processor.executor >> executor
        executor.batchService >> batchService
        batchService.config >> config
        config.batch() >> batchOpts
        batchOpts.terminateJobsOnCompletion >> false

        when:
        def observer = new AzBatchProcessObserver(Mock(Session))
        observer.onProcessTerminate(processor)

        then:
        0 * batchService.getJobIdsForProcessor(_)
        0 * batchService.setJobTermination(_)
    }

    def 'should handle exception from setJobTermination gracefully'() {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(AzBatchExecutor)
        def batchService = Mock(AzBatchService)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts)
        and:
        processor.executor >> executor
        processor.name >> 'test_process'
        executor.batchService >> batchService
        batchService.config >> config
        config.batch() >> batchOpts
        batchOpts.terminateJobsOnCompletion >> true
        batchService.getJobIdsForProcessor(processor) >> ['job-1', 'job-2']
        batchService.setJobTermination('job-1') >> { throw new RuntimeException('API error') }

        when:
        def observer = new AzBatchProcessObserver(Mock(Session))
        observer.onProcessTerminate(processor)

        then:
        noExceptionThrown()
        1 * batchService.setJobTermination('job-2')
    }

    def 'should handle empty job list'() {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(AzBatchExecutor)
        def batchService = Mock(AzBatchService)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts)
        and:
        processor.executor >> executor
        processor.name >> 'test_process'
        executor.batchService >> batchService
        batchService.config >> config
        config.batch() >> batchOpts
        batchOpts.terminateJobsOnCompletion >> true
        batchService.getJobIdsForProcessor(processor) >> []

        when:
        def observer = new AzBatchProcessObserver(Mock(Session))
        observer.onProcessTerminate(processor)

        then:
        0 * batchService.setJobTermination(_)
    }
}
