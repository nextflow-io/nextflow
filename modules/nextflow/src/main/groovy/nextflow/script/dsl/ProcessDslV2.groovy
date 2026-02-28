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

package nextflow.script.dsl

import groovy.transform.TypeChecked
import nextflow.script.BaseScript
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.params.v2.ProcessFileInput
import nextflow.script.params.v2.ProcessFileOutput
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessInputsDef
import nextflow.script.params.v2.ProcessOutputsDef

/**
 * Implements the DSL for typed processes.
 *
 * @see nextflow.script.control.ScriptToGroovyVisitor
 * @see nextflow.script.dsl.ProcessDsl.StageDsl
 * @see nextflow.script.dsl.ProcessDsl.OutputDslV2
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@TypeChecked
class ProcessDslV2 extends ProcessBuilder {

    private ProcessInputsDef inputs

    private ProcessOutputsDef outputs

    ProcessDslV2(BaseScript ownerScript, String processName) {
        super(new ProcessConfigV2(ownerScript, processName))
        inputs = ((ProcessConfigV2) config).getInputs()
        outputs = ((ProcessConfigV2) config).getOutputs()
    }

    /// INPUTS

    /**
     * Declare a process input.
     *
     * @param name
     * @param type
     * @param optional
     */
    void _input_(String name, Class type, boolean optional) {
        inputs.addParam(name, type, optional)
    }

    /**
     * Declare a process tuple input.
     *
     * @param components
     * @param type
     */
    void _input_(List<ProcessInput> components, Class type) {
        inputs.addTupleParam(components, type)
    }

    /// STAGE DIRECTIVES

    /**
     * Declare an environment variable in the task environment
     * with the given name and value.
     *
     * @paran name
     * @param value [String | Closure]
     */
    void env(String name, Object value) {
        inputs.addEnv(name, value)
    }

    /**
     * Declare a file or collection of files to be staged into
     * the task directory.
     *
     * This method is automatically generated for Path inputs.
     *
     * @param value [Path | Collection<Path> | Closure]
     */
    void stageAs(Object value) {
        inputs.addFile(new ProcessFileInput(null, value))
    }

    /**
     * Declare a file or collection of files to be staged into
     * the task directory under the given file pattern.
     *
     * @param filePattern [String | Closure]
     * @param value       [Path | Collection<Path> | Closure]
     */
    void stageAs(Object filePattern, Object value) {
        inputs.addFile(new ProcessFileInput(filePattern, value))
    }

    /**
     * Declare a value as the standard input of the task script.
     *
     * @param value [String | Closure]
     */
    void stdin(Object value) {
        inputs.setStdin(value)
    }

    /// OUTPUTS

    /**
     * Declare a process output.
     *
     * @param name
     * @param type
     * @param value
     */
    void _output_(String name, Class type, Object value) {
        outputs.addParam(name, type, value)
    }

    /**
     * Declare a value to be emitted to the given topic
     * on task completion.
     *
     * @param value
     * @param target
     */
    void _topic_(Object value, String target) {
        outputs.addTopic(value, target)
    }

    /// UNSTAGE DIRECTIVES

    /**
     * Declare an output environment variable to be extracted from
     * the task environment.
     *
     * @param name
     */
    void _unstage_env(String name) {
        outputs.addEnv(name)
    }

    /**
     * Declare an eval command to be executed in the task environment.
     *
     * @param key
     * @param cmd [String | Closure<String>]
     */
    void _unstage_eval(String key, Object cmd) {
        outputs.addEval(key, cmd)
    }

    /**
     * Declare an output file or collection of files to be unstaged
     * from the task environment.
     *
     * @param key
     * @param pattern [String | Closure<String>]
     */
    void _unstage_files(String key, Object pattern) {
        outputs.addFile(key, new ProcessFileOutput(pattern))
    }

}
