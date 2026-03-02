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

package nextflow.script

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
/**
 * Models Seqera Platform metadata for Nextflow execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
class PlatformMetadata {

    @Canonical
    static class User {
        String id
        String userName
        String email
        String firstName
        String lastName
        String organization

        User(Map opts) {
            id = opts.id as String
            userName = opts.userName as String
            email = opts.email as String
            firstName = opts.firstName as String
            lastName = opts.lastName as String
            organization = opts.organization as String
        }
    }

    @Canonical
    static class Workspace {
        String id
        String name
        String fullName
        String organization

        Workspace(Map opts) {
            id = opts.workspaceId as String
            name = opts.workspaceName as String
            fullName = opts.workspaceFullName
            organization = opts.orgName as String
        }
    }

    @Canonical
    static class ComputeEnv {
        String id
        String name
        String platform

        ComputeEnv(Map opts) {
            id = opts.id as String
            name = opts.name as String
            platform = opts.platform as String
        }
    }

    @Canonical
    static class Pipeline {
        String id
        String name
        String revision
        String commitId

        Pipeline(Map opts) {
            id = opts.id as String
            name = opts.name as String
            revision = opts.revision as String
            commitId = opts.commitId as String
        }
    }

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
    
    volatile User user
    volatile Workspace workspace
    volatile ComputeEnv computeEnv
    volatile Pipeline pipeline
    volatile List labels

    PlatformMetadata() {}

    PlatformMetadata(String workflowId) {
        this.workflowId = workflowId
    }
    
}
