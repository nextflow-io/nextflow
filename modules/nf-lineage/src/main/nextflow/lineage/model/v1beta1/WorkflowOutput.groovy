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

package nextflow.lineage.model.v1beta1

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.lineage.serde.LinSerializable

import java.time.OffsetDateTime

/**
 * Models the results of a workflow execution.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io
 */
@Canonical
@CompileStatic
class WorkflowOutput implements LinSerializable {
    /**
     * Creation date of the workflow output
     */
    OffsetDateTime createdAt
    /**
     * Workflow run that generated the output
     */
    String workflowRun
    /**
     * Workflow output
     */
    List<Parameter> output
}
