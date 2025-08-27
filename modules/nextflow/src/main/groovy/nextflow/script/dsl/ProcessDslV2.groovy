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

package nextflow.script.dsl

import groovy.transform.CompileStatic
import nextflow.script.BaseScript
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.params.v2.ProcessFileInput
import nextflow.script.params.v2.ProcessFileOutput
import nextflow.script.params.v2.ProcessInputs
import nextflow.script.params.v2.ProcessOutputs

/**
 * Implements the DSL for typed processes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessDslV2 extends ProcessBuilder {

    private ProcessInputs inputs = new ProcessInputs()

    private ProcessOutputs outputs = new ProcessOutputs()

    ProcessDslV2(BaseScript ownerScript, String processName) {
        super(new ProcessConfigV2(ownerScript, processName))
    }

    /// INPUTS

    void _input_(String name, Class type) {
        inputs.addParam(name, type)
    }

    /// STAGE DIRECTIVES

    void env(String name, Object value) {
        inputs.addEnv(name, value)
    }

    void stageAs(Object value) {
        inputs.addFile(new ProcessFileInput(null, value))
    }

    void stageAs(Object name, Object value) {
        inputs.addFile(new ProcessFileInput(name, value))
    }

    void stdin(Object value) {
        inputs.setStdin(value)
    }

    /// OUTPUTS

    void _output_(String name, Object value) {
        outputs.addParam(name, null, value)
    }

    void _topic_(Object value, String topic) {
        // TODO
    }

    /// UNSTAGE DIRECTIVES

    void _unstage_env(String name) {
        outputs.addEnv(name)
    }

    void _unstage_eval(String key, String cmd) {
        outputs.addEval(key, cmd)
    }

    void _unstage_files(String key, Object pattern) {
        outputs.addFile(key, new ProcessFileOutput(pattern))
    }

    /// BUILD

    ProcessDef build() {
        ((ProcessConfigV2) config).setInputs(inputs)
        ((ProcessConfigV2) config).setOutputs(outputs)
        super.build()
    }

}
