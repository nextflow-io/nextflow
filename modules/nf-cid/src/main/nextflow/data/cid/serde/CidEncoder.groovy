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

package nextflow.data.cid.serde


import groovy.transform.CompileStatic
import nextflow.data.cid.model.TaskOutput
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutput
import nextflow.data.cid.model.WorkflowResults
import nextflow.data.cid.model.WorkflowRun
import nextflow.serde.gson.GsonEncoder
import nextflow.serde.gson.RuntimeTypeAdapterFactory
/**
 * Implements a JSON encoder for CID model objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidEncoder extends GsonEncoder<CidSerializable> {
    public static RuntimeTypeAdapterFactory CID_SERIALIZABLE_FACTORY = RuntimeTypeAdapterFactory.of(CidSerializable.class, "type")
        .registerSubtype(WorkflowRun, WorkflowRun.simpleName)
        .registerSubtype(WorkflowResults, WorkflowResults.simpleName)
        .registerSubtype(Workflow, Workflow.simpleName)
        .registerSubtype(WorkflowOutput, WorkflowOutput.simpleName)
        .registerSubtype(TaskRun, TaskRun.simpleName)
        .registerSubtype(TaskOutput, TaskOutput.simpleName)

    CidEncoder() {
        withTypeAdapterFactory(CID_SERIALIZABLE_FACTORY)
    }

}
