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
import nextflow.script.ProcessFileOutput
import nextflow.script.ProcessOutput
import nextflow.script.ProcessOutputs

/**
 * Builder for {@link ProcessOutputs}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessOutputsBuilder {

    private ProcessOutputs outputs = new ProcessOutputs()

    ProcessOutputsBuilder env(String name) {
        env(name, name)
    }

    ProcessOutputsBuilder env(String name, String target) {
        outputs.env.put(name, target)
        return this
    }

    ProcessOutputsBuilder path(Map opts=[:], String name) {
        path(opts, name, name)
    }

    ProcessOutputsBuilder path(Map opts=[:], String name, Object target) {
        outputs.files.put(name, new ProcessFileOutput(target, opts))
        return this
    }

    /**
     * Declare a process output with a closure that will
     * be evaluated after the task execution. For example:
     *
     *   env 'SAMPLE_ID'           // declare output env 'SAMPLE_ID'
     *   path '$out0', 'file.txt'  // declare output file 'file.txt'
     *
     *   emit { sample_id }        // variable 'sample_id' in task context
     *   emit { stdout }           // standard output of task script
     *   emit { [env('SAMPLE_ID'), path('$out0')] }
     *   emit { new Sample(sample_id, path('$out0')) }
     *
     * @param opts
     * @param target
     */
    ProcessOutputsBuilder emit(Map opts=[:], Object target) {
        outputs.add(new ProcessOutput(outputs, target, opts))
        return this
    }

    ProcessOutputs build() {
        outputs
    }

}
