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
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.script.params.InParam

/**
 * Models a process input declaration.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessInput implements InParam {

    /**
     * Parameter name under which the input value for each task
     * will be added to the task context.
     */
    private String name

    /**
     * Parameter type which is used to validate task inputs
     */
    private Class type

    /**
     * Whether the input can be null.
     */
    private boolean optional

    /**
     * Input channel which is created when the process is invoked
     * in a workflow.
     */
    DataflowReadChannel channel

    ProcessInput(String name, Class type, boolean optional) {
        this.name = name
        this.type = type
        this.optional = optional
    }

    @Override
    String getName() {
        return name
    }

    Class getType() {
        return type
    }

    boolean isOptional() {
        return optional
    }

    @Override
    ProcessInput clone() {
        (ProcessInput)super.clone()
    }

    /// LEGACY METHODS

    @Override
    DataflowReadChannel getInChannel() {
        return channel
    }

    @Override
    Object getRawChannel() {
        throw new UnsupportedOperationException()
    }

    @Override
    def decodeInputs( List values ) {
        throw new UnsupportedOperationException()
    }

}
