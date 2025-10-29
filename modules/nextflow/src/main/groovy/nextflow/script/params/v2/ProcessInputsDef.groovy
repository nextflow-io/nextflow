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
import nextflow.extension.CH

/**
 * Models the process `input` section, including the input
 * declarations and staging directives.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessInputsDef implements Cloneable {

    private List<ProcessInput> params = []

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

    void addParam(String name, Class type, boolean optional) {
        params.add(new ProcessInput(name, type, optional))
    }

    void addTupleParam(List<ProcessInput> components, Class type) {
        params.add(new ProcessTupleInput(components, type))
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

    ProcessInput getAt(int i) {
        return params.get(i)
    }

    Map<String,?> getEnv() {
        return env
    }

    List<ProcessFileInput> getFiles() {
        return files
    }

    List<DataflowReadChannel> getChannels() {
        final result = new ArrayList<DataflowReadChannel>()
        for( final param : params )
            result.add(param.getChannel())
        return result
    }

    boolean isSingleton() {
        return params.every { param ->
            !CH.isChannelQueue(param.getChannel())
        }
    }

    @Override
    ProcessInputsDef clone() {
        final result = (ProcessInputsDef)super.clone()
        result.params = new ArrayList<>(params.size())
        for( final param : params ) {
            result.params.add(param.clone())
        }
        return result
    }

}