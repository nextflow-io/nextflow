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

package nextflow.script

import nextflow.script.params.DefaultInParam
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutputsList

/**
 * Specialization of process config for legacy processes.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessConfigV1 extends ProcessConfig {

    private InputsList inputs = new InputsList()

    private OutputsList outputs = new OutputsList()

    ProcessConfigV1(BaseScript script, String name) {
        super(script, name)
    }

    @Override
    ProcessConfigV1 clone() {
        final copy = (ProcessConfigV1)super.clone()
        copy.@inputs = inputs.clone()
        copy.@outputs = outputs.clone()
        return copy
    }

    InputsList getInputs() {
        return inputs
    }

    OutputsList getOutputs() {
        return outputs
    }

    /**
     * Defines a special *dummy* input parameter, when no inputs are
     * provided by the user for the current task
     */
    void fakeInput() {
        new DefaultInParam(this)
    }

    void fakeOutput() {
        new DefaultOutParam(this)
    }

}
