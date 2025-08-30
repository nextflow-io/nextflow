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

package nextflow.script

import nextflow.script.params.v2.ProcessInputs
import nextflow.script.params.v2.ProcessOutputs

/**
 * Specialization of process config for typed processes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessConfigV2 extends ProcessConfig {

    /**
     * List of process input definitions
     */
    private ProcessInputs inputs = new ProcessInputs()

    /**
     * List of process output definitions
     */
    private ProcessOutputs outputs = new ProcessOutputs()

    ProcessConfigV2(BaseScript script, String name) {
        super(script, name)
    }

    @Override
    ProcessConfigV2 clone() {
        final copy = (ProcessConfigV2)super.clone()
        copy.@inputs = inputs.clone()
        copy.@outputs = outputs.clone()
        return copy
    }

    ProcessInputs getInputs() {
        return inputs
    }

    ProcessOutputs getOutputs() {
        return outputs
    }

}
