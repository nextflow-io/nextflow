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

    def 'should only act on Azure Batch executors'() {
        given:
        def observer = new AzBatchProcessObserver(Mock(Session))
        def processor = Mock(TaskProcessor) {
            getExecutor() >> Mock(Executor)  // Not an AzBatchExecutor
        }

        when:
        observer.onProcessTerminate(processor)

        then:
        noExceptionThrown()
    }

    def 'should set job auto-termination when enabled'() {
        given:
        def observer = new AzBatchProcessObserver(Mock(Session))
        def processor = Mock(TaskProcessor) { getName() >> 'test-process' }
        def batchService = Mock(AzBatchService) {
            getConfig() >> Mock(AzConfig) {
                batch() >> Mock(AzBatchOpts) {
                    terminateJobsOnCompletion >> true
                }
            }
            getAllJobIds() >> [(new AzJobKey(processor, 'pool1')): 'job123']
        }
        def executor = Mock(AzBatchExecutor) { getBatchService() >> batchService }
        processor.getExecutor() >> executor

        when:
        observer.onProcessTerminate(processor)

        then:
        1 * batchService.setJobAutoTermination('job123')
    }

    def 'should skip when termination disabled'() {
        given:
        def observer = new AzBatchProcessObserver(Mock(Session))
        def processor = Mock(TaskProcessor) { getName() >> 'test-process' }
        def batchService = Mock(AzBatchService) {
            getConfig() >> Mock(AzConfig) {
                batch() >> Mock(AzBatchOpts) {
                    terminateJobsOnCompletion >> false
                }
            }
        }
        def executor = Mock(AzBatchExecutor) { getBatchService() >> batchService }
        processor.getExecutor() >> executor

        when:
        observer.onProcessTerminate(processor)

        then:
        0 * batchService.setJobAutoTermination(_)
    }
} 