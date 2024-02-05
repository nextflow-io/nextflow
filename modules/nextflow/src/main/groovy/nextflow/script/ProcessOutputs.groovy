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

package nextflow.script

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.util.LazyVar

/**
 * Models the process outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessOutputs implements List<ProcessOutput>, Cloneable {

    @Delegate
    private List<ProcessOutput> params = []

    /**
     * Environment variables which will be exported from the
     * task environment for each task and made available to
     * process outputs.
     */
    private Map<String,String> env = [:]

    /**
     * Output files which will be unstaged from the task
     * directory for each task and made available to process
     * outputs.
     */
    private Map<String,ProcessFileOutput> files = [:]

    @Override
    ProcessOutputs clone() {
        def result = (ProcessOutputs)super.clone()
        result.params = new ArrayList<>(params.size())
        for( ProcessOutput param : params )
            result.add((ProcessOutput)param.clone())
        return result
    }

    void addParam(Object target, Map opts) {
        add(new ProcessOutput(this, target, opts))
    }

    void setDefault() {
        final param = new ProcessOutput(this, new LazyVar('stdout'), [:])
        param.setChannel(new DataflowQueue())
        params.add(param)
    }

    void addEnv(String name, Object value) {
        env.put(name, value)
    }

    void addFile(String key, ProcessFileOutput file) {
        files.put(key, file)
    }

    List<String> getNames() {
        return params*.getName()
    }

    List<DataflowWriteChannel> getChannels() {
        return params*.getChannel()
    }

    Map<String,String> getEnv() {
        return env
    }

    Map<String,ProcessFileOutput> getFiles() {
        return files
    }

}
