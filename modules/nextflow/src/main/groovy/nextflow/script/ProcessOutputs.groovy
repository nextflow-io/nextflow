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

/**
 * Models the process outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessOutputs implements List<ProcessOutput>, Cloneable {

    @Delegate
    private List<ProcessOutput> target = []

    Map<String,String> env = [:]

    Map<String,ProcessFileOutput> files = [:]

    @Override
    ProcessOutputs clone() {
        def result = (ProcessOutputs)super.clone()
        result.target = new ArrayList<>(target.size())
        for( ProcessOutput param : target )
            result.add((ProcessOutput)param.clone())
        return result
    }

    List<String> getNames() {
        return target*.getName()
    }

    List<DataflowWriteChannel> getChannels() {
        return target*.getChannel()
    }

    void setDefault() {
        final param = new ProcessOutput(ProcessOutput.Shortcuts.STDOUT, [:])
        param.setChannel(new DataflowQueue())
        target.add(param)
    }

}
