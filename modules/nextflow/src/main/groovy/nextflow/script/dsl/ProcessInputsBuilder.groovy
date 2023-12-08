/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.script.dsl

import groovy.transform.CompileStatic
import nextflow.script.ProcessFileInput
import nextflow.script.ProcessInput
import nextflow.script.ProcessInputs

/**
 * Builder for {@link ProcessInputs}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessInputsBuilder {

    private ProcessInputs inputs = new ProcessInputs()

    ProcessInputsBuilder env(String name, Object source) {
        inputs.env.put(name, source)
        return this
    }

    ProcessInputsBuilder path(Map opts=[:], Object source) {
        inputs.files.add(new ProcessFileInput(source, null, true, opts))
        return this
    }

    ProcessInputsBuilder stdin(Object source) {
        inputs.stdin = source
        return this
    }

    /**
     * Declare a process input parameter which will be
     * bound when the task is created and can be referenced by
     * other process input methods. For example:
     *
     *   take 'sample_id'
     *   take 'files'
     *
     *   env('SAMPLE_ID') { sample_id }
     *   path { files }
     *
     * @param name
     */
    ProcessInputsBuilder take(String name) {
        inputs.add(new ProcessInput(name))
        return this
    }

    ProcessInputs build() {
        return inputs
    }

}
