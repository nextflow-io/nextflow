/*
 * Copyright 2021, Microsoft Corp
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
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Azure extension of the TaskBean that includes Azure Batch pool options
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzTaskBean extends TaskBean {
    
    /**
     * The Azure pool options that apply to this task
     */
    final AzPoolOpts poolOpts
    
    /**
     * Create a new Azure TaskBean from a TaskRun with pool options
     * 
     * @param task The TaskRun to get base configuration from
     * @param poolOpts The Azure pool options
     */
    AzTaskBean(TaskRun task, AzPoolOpts poolOpts) {
        super(task)
        this.poolOpts = poolOpts
    }
    
    /**
     * Get the Azure pool options
     * 
     * @return The Azure pool options
     */
    AzPoolOpts getPoolOpts() {
        return poolOpts
    }
    
    /**
     * Get the managed identity client ID if configured
     * 
     * @return The managed identity client ID or null if not configured
     */
    String getManagedIdentityId() {
        return poolOpts?.managedIdentityId
    }
} 