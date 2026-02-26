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

package nextflow.script.params.v2

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.script.params.OutParam

/**
 * Models a process output declaration.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessOutput implements OutParam {

    /**
     * Name of the output channel in the process outputs (i.e. `.out`).
     */
    private String name

    /**
     * Optional output type which is used to validate task outputs.
     */
    private Class type

    /**
     * Lazy expression (e.g. closure) which defines the output value
     * in terms of the task context, including environment variables,
     * files, and standard output.
     * It will be evaluated for each task after it is executed. 
     */
    private Object value

    /**
     * Output channel which is created when the process is invoked
     * in a workflow.
     */
    DataflowWriteChannel channel

    ProcessOutput(String name, Class type, Object value) {
        this.name = name
        this.type = type
        this.value = value
    }

    @Override
    String getName() {
        return name
    }

    Object getLazyValue() {
        return value
    }

    @Override
    ProcessOutput clone() {
        (ProcessOutput)super.clone()
    }

    /// LEGACY METHODS

    DataflowWriteChannel getOutChannel() {
        return channel
    }

    short getIndex() {
        throw new UnsupportedOperationException()
    }

    String getChannelEmitName() {
        throw new UnsupportedOperationException()
    }

    String getChannelTopicName() {
        throw new UnsupportedOperationException()
    }

}
