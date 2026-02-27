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

package nextflow.lineage.model.v1beta1

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.lineage.serde.LinSerializable

/**
 * Models a task execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class TaskRun implements LinSerializable {
    /**
     * Execution session identifier
     */
    String sessionId
    /**
     * Task name
     */
    String name
    /**
     * Checksum of the task source code
     */
    Checksum codeChecksum
    /**
     * Resolved task script
     */
    String script
    /**
     * Task run input
     */
    List<Parameter> input
    /**
     * Container used for the task run
     */
    String container
    /**
     * Conda environment used for the task run
     */
    String conda
    /**
     * Spack environment used for the task run
     */
    String spack
    /**
     * Architecture defined in the Spack environment used for the task run
     */
    String architecture
    /**
     * Global variables defined in the task run
     */
    Map globalVars
    /**
     * Binaries used in the task run
     */
    List<DataPath> binEntries
    /**
     * Workflow run associated to the task run
     */
    String workflowRun
}
