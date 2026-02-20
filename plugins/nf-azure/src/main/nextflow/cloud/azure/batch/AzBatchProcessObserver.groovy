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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.trace.TraceObserver

/**
 * Observer that handles process termination events for Azure Batch executor.
 * When a process terminates (all tasks have been submitted and completed), this observer
 * will set the corresponding Azure Batch jobs to auto-terminate when all tasks complete.
 *
 * This provides an additional mechanism to ensure jobs are terminated during workflow
 * execution, complementing the eager per-task termination and the workflow-end fallback.
 *
 * @author Adam Talbot <adam.talbot@seqera.io>
 */
@Slf4j
@CompileStatic
class AzBatchProcessObserver implements TraceObserver {

    AzBatchProcessObserver(Session session) {
        // Session not needed, but required by factory interface
    }

    /**
     * Called when a process terminates (all tasks have been submitted and completed).
     * Sets Azure Batch jobs to auto-terminate when all tasks complete.
     */
    @Override
    void onProcessTerminate(TaskProcessor processor) {
        // Check if this process uses the Azure Batch executor
        if( !(processor.executor instanceof AzBatchExecutor) ) {
            return
        }

        final executor = processor.executor as AzBatchExecutor
        final batchService = executor.batchService

        // Check if auto-termination is enabled
        if( !batchService?.config?.batch()?.terminateJobsOnCompletion ) {
            log.trace "Azure Batch job auto-termination is disabled, skipping process termination for: ${processor.name}"
            return
        }

        // Find and set auto-termination for all jobs associated with this processor
        batchService.getJobIdsForProcessor(processor).each { jobId ->
            log.debug "Setting Azure Batch job ${jobId} to auto-terminate for completed process: ${processor.name}"
            try {
                batchService.setJobTermination(jobId)
            }
            catch( Exception e ) {
                log.debug "Failed to set auto-termination for Azure Batch job ${jobId} associated with process '${processor.name}' - ${e.message ?: e}"
            }
        }
    }
}
