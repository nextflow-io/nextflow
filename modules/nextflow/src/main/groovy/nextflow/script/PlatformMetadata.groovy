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

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Models Seqera Platform metadata for Nextflow execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
class PlatformMetadata {

    /**
     * Volatile because it is written by TowerClient.onFlowCreate on the main thread
     * and read by SeqeraExecutor.createRun on the executor thread.
     */
    volatile String workflowId

    /**
     * The Platform watch URL for the current workflow execution.
     * Set by TowerClient.onFlowBegin and read by SeqeraExecutor.createRun.
     */
    volatile String workflowUrl

    PlatformMetadata() {}

    PlatformMetadata(String workflowId) {
        this.workflowId = workflowId
    }
}
