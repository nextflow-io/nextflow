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
 */

package nextflow.script

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j

/**
 * Store the workflow platform-related metadata
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
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

    String workflowId
    User user
    Workspace workspace
    ComputeEnv computeEnv
    Pipeline pipeline
    List labels

    PlatformMetadata() {}

    PlatformMetadata(String id, User user, Workspace workspace, ComputeEnv computeEnv, Pipeline pipeline, List labels) {
        this.workflowId = id
        this.user = user
        this.workspace = workspace
        this.computeEnv = computeEnv
        this.pipeline = pipeline
        this.labels = labels
    }
}