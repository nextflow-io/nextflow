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
import nextflow.script.params.InParam

/**
 * Models the process inputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessInputs implements Cloneable {

    private List<ProcessInput> params = []

    /**
     * Input variables which will be evaluated for each task
     * in terms of the task inputs and added to the task context.
     */
    private Map<String,?> vars = [:]

    /**
     * Environment variables which will be evaluated for each
     * task against the task context and added to the task
     * environment.
     */
    private Map<String,?> env = [:]

    /**
     * Input files which will be evaluated for each task
     * against the task context and staged into the task
     * directory.
     */
    private List<ProcessFileInput> files = []

    /**
     * Lazy expression which will be evaluated for each task
     * against the task context and provided as the standard
     * input to the task.
     */
    Object stdin

    /**
     * Input channel which is created when the process is invoked
     * in a workflow.
     */
    DataflowReadChannel channel

    void addParam(String name, Class type) {
        params.add(new ProcessInput(name, type))
    }

    void addVariable(String name, Object value) {
        vars.put(name, value)
    }

    void addEnv(String name, Object value) {
        env.put(name, value)
    }

    void addFile(ProcessFileInput file) {
        files.add(file)
    }

    List<ProcessInput> getParams() {
        return params
    }

    int size() {
        return params.size()
    }

    Map<String,?> getVariables() {
        return vars
    }

    Map<String,?> getEnv() {
        return env
    }

    List<ProcessFileInput> getFiles() {
        return files
    }

    @Override
    ProcessInputs clone() {
        final result = (ProcessInputs)super.clone()
        result.params = new ArrayList<>(params.size())
        for( final param : params ) {
            result.params.add(param.clone())
        }
        return result
    }

}


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

    ProcessInput(String name, Class type) {
        this.name = name
        this.type = type
    }

    @Override
    String getName() {
        return name
    }

    Class getType() {
        return type
    }

    @Override
    ProcessInput clone() {
        (ProcessInput)super.clone()
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
