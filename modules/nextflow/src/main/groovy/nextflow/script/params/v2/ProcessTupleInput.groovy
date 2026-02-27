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
import groovyx.gpars.dataflow.DataflowReadChannel

/**
 * Models a process tuple input, which declares a separate
 * input for each tuple component.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessTupleInput extends ProcessInput {

    private List<ProcessInput> components

    ProcessTupleInput(List<ProcessInput> components, Class type) {
        super("", type, false)
        this.components = components
    }

    List<ProcessInput> getComponents() {
        return components
    }

    @Override
    ProcessTupleInput clone() {
        (ProcessTupleInput)super.clone()
    }

    /// LEGACY METHODS

    @Override
    DataflowReadChannel getInChannel() {
        throw new UnsupportedOperationException()
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
