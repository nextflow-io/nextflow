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

package nextflow.data.cid.model

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.data.cid.serde.CidSerializable

import java.time.Instant

/**
 * Models the results of a workflow execution.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io
 */
@Canonical
@CompileStatic
class WorkflowOutputs implements CidSerializable {
    /**
     * Creation date of the workflow outputs description
     */
    Instant createdAt
    /**
     * Workflow run that generated the outputs
     */
    String workflowRun
    /**
     * Workflow outputs
     */
    Map<String, Object> outputs
    /**
     * Annotations attached to the workflow outputs
     */
    Map annotations
}
