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
 * When a process terminates (all tasks have been submitted), this observer 
 * will eagerly set the corresponding Azure Batch job to auto-terminate when
 * all tasks complete.
 *
 * @author Adam Talbot <adam.talbot@seqera.io>
 */
@Slf4j
@CompileStatic
class AzBatchProcessObserver implements TraceObserver {

    private Session session

    AzBatchProcessObserver(Session session) {
        this.session = session
    }

    /**
     * Called when a process terminates (all tasks have been submitted).
     * This method will find the Azure Batch job ID associated with the process
     * and set it to auto-terminate when all tasks complete.
     */
    @Override
    void onProcessTerminate(TaskProcessor processor) {
        // Check if this process uses the Azure Batch executor
        if( !(processor.executor instanceof AzBatchExecutor) ) {
            return
        }
        
        final executor = processor.executor as AzBatchExecutor
        final batchService = executor.getBatchService()
        if( !batchService?.config?.batch()?.terminateJobsOnCompletion ) {
            log.trace "Azure Batch job auto-termination is disabled, skipping eager termination for process: ${processor.name}"
            return
        }

        try {
            // Find all job IDs associated with this processor
            final jobIds = batchService.allJobIds.findAll { key, jobId ->
                key.processor == processor
            }.values()

            for( String jobId : jobIds ) {
                log.debug "Setting Azure Batch job ${jobId} to auto-terminate for completed process: ${processor.name}"
                batchService.setJobAutoTermination(jobId)
            }
        }
        catch( Exception e ) {
            log.warn "Failed to set auto-termination for Azure Batch jobs associated with process '${processor.name}' - Reason: ${e.message ?: e}"
        }
    }
} 