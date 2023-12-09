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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.TaskOutputCollector
import nextflow.processor.TaskRun
import nextflow.util.LazyHelper
/**
 * Models a process output.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessOutput implements Cloneable {

    private ProcessOutputs declaredOutputs

    private Object target

    private String name

    private boolean optional

    private DataflowWriteChannel channel

    ProcessOutput(ProcessOutputs declaredOutputs, Object target, Map<String,?> opts) {
        this.declaredOutputs = declaredOutputs
        this.target = target

        if( opts.name )
            this.name = opts.name
        if( opts.optional )
            this.optional = true
    }

    String getName() {
        return name
    }

    void setChannel(DataflowWriteChannel channel) {
        this.channel = channel
    }

    DataflowWriteChannel getChannel() {
        return channel
    }

    Object resolve(TaskRun task) {
        final ctx = new TaskOutputCollector(declaredOutputs, optional, task)
        return LazyHelper.resolve(ctx, target)
    }

}
