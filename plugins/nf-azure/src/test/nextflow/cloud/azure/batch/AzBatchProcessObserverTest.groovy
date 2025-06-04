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
import nextflow.cloud.azure.config.AzBatchOpts
import nextflow.cloud.azure.config.AzConfig
import nextflow.executor.Executor
import nextflow.processor.TaskProcessor
import spock.lang.Specification

/**
 * Test for AzBatchProcessObserver
 *
 * @author Adam Talbot <adam.talbot@seqera.io>
 */
class AzBatchProcessObserverTest extends Specification {

    def 'should skip non-azure batch executors'() {
        given:
        def session = Mock(Session)
        def observer = new AzBatchProcessObserver(session)
        def executor = Mock(Executor)  // Not an AzBatchExecutor
        def processor = Mock(TaskProcessor) {
            getExecutor() >> executor
        }

        when:
        observer.onProcessTerminate(processor)

        then:
        // No exception should be thrown and method should return early
        noExceptionThrown()
    }

    def 'should skip when termination is disabled'() {
        given:
        def session = Mock(Session)
        def observer = new AzBatchProcessObserver(session)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts) {
            terminateJobsOnCompletion >> false
        }
        def batchService = Mock(AzBatchService) {
            getConfig() >> config
        }
        def executor = Mock(AzBatchExecutor) {
            getBatchService() >> batchService
        }
        def processor = Mock(TaskProcessor) {
            getExecutor() >> executor
            getName() >> 'test-process'
        }
        config.batch() >> batchOpts

        when:
        observer.onProcessTerminate(processor)

        then:
        0 * batchService.setJobAutoTermination(_)
    }

    def 'should set job auto-termination when enabled'() {
        given:
        def session = Mock(Session)
        def observer = new AzBatchProcessObserver(session)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts) {
            terminateJobsOnCompletion >> true
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'test-process'
        }
        def jobKey1 = new AzJobKey(Mock(TaskProcessor), 'pool1')
        def jobKey2 = new AzJobKey(processor, 'pool1')  // This one should match
        def jobKey3 = new AzJobKey(Mock(TaskProcessor), 'pool2')
        def allJobIds = [
            (jobKey1): 'job1',
            (jobKey2): 'job2',  // This should be processed
            (jobKey3): 'job3'
        ]
        def batchService = Mock(AzBatchService) {
            getConfig() >> config
            getAllJobIds() >> allJobIds
        }
        def executor = Mock(AzBatchExecutor) {
            getBatchService() >> batchService
        }
        processor.getExecutor() >> executor
        config.batch() >> batchOpts

        when:
        observer.onProcessTerminate(processor)

        then:
        1 * batchService.setJobAutoTermination('job2')
    }

    def 'should handle exceptions gracefully'() {
        given:
        def session = Mock(Session)
        def observer = new AzBatchProcessObserver(session)
        def config = Mock(AzConfig)
        def batchOpts = Mock(AzBatchOpts) {
            terminateJobsOnCompletion >> true
        }
        def batchService = Mock(AzBatchService) {
            getConfig() >> config
            getAllJobIds() >> { throw new RuntimeException('Test error') }
        }
        def executor = Mock(AzBatchExecutor) {
            getBatchService() >> batchService
        }
        def processor = Mock(TaskProcessor) {
            getExecutor() >> executor
            getName() >> 'test-process'
        }
        config.batch() >> batchOpts

        when:
        observer.onProcessTerminate(processor)

        then:
        // Should not throw exception
        noExceptionThrown()
    }
} 