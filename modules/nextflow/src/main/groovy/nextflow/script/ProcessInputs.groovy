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

import groovyx.gpars.dataflow.DataflowReadChannel

/**
 * Models the process inputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessInputs implements List<ProcessInput>, Cloneable {

    @Delegate
    private List<ProcessInput> target = []

    Map<String,?> env = [:]

    List<ProcessFileInput> files = []

    Object stdin

    @Override
    ProcessInputs clone() {
        def result = (ProcessInputs)super.clone()
        result.target = new ArrayList<>(target.size())
        for( ProcessInput param : target ) {
            result.target.add((ProcessInput)param.clone())
        }
        return result
    }

    List<String> getNames() {
        return target*.getName()
    }

    List<DataflowReadChannel> getChannels() {
        return target*.getChannel()
    }

}
