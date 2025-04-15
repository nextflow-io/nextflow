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
import java.time.OffsetDateTime

/**
 * Models task results.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Canonical
@CompileStatic
class TaskOutputs implements CidSerializable {
    /**
     * Reference to the task that generated the data.
     */
    String taskRun
    /**
     * Reference to the WorkflowRun that generated the data.
     */
    String workflowRun
    /**
     * Creation date of this task outputs description
     */
    OffsetDateTime createdAt
    /**
     * Outputs of the task
     */
    List<Parameter> outputs
    /**
     * Annotations attached to the task outputs
     */
    List<Annotation> annotations
}
