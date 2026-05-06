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

package nextflow.script.params.v2

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.script.params.OutParam

/**
 * Models a process topic emission.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessTopic {

    /**
     * Lazy expression (e.g. closure) which defines the output value
     * in terms of the task context, including environment variables,
     * files, and standard output.
     * It will be evaluated for each task after it is executed.
     */
    private Object value

    /**
     * Name of the target topic.
     */
    private String target

    /**
     * Topic channel which is bound when the process is invoked
     * in a workflow.
     */
    DataflowWriteChannel channel

    ProcessTopic(Object value, String target) {
        this.value = value
        this.target = target
    }

    Object getLazyValue() {
        return value
    }

    String getTarget() {
        return target
    }

    @Override
    ProcessTopic clone() {
        (ProcessTopic)super.clone()
    }

}
