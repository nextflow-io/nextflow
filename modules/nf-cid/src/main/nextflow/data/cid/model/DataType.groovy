/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.data.cid.model

/**
 * Possible metadata type entries.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
enum DataType {
    TaskRun(nextflow.data.cid.model.TaskRun),
    Workflow(nextflow.data.cid.model.Workflow),
    WorkflowRun(nextflow.data.cid.model.WorkflowRun),
    TaskOutput(nextflow.data.cid.model.Output),
    WorkflowOutput(nextflow.data.cid.model.Output),
    WorkflowResults(nextflow.data.cid.model.WorkflowResults)

    final Class clazz

    DataType(Class clazz) {
        this.clazz = clazz
    }
}
