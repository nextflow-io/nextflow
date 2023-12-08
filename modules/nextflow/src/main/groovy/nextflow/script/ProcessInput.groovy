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
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.extension.CH

/**
 * Models a process input.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessInput implements Cloneable {

    private String name

    private Object arg

    private DataflowReadChannel channel

    ProcessInput(String name) {
        this.name = name
    }

    String getName() {
        return name
    }

    void bind(Object arg) {
        this.arg = arg
        this.channel = getInChannel(arg)
    }

    private DataflowReadChannel getInChannel(Object obj) {
        if( obj == null )
            throw new IllegalArgumentException('A process input channel evaluates to null')

        final result = obj instanceof Closure
            ? obj.call()
            : obj

        if( result == null )
            throw new IllegalArgumentException('A process input channel evaluates to null')

        def inChannel
        if ( result instanceof DataflowReadChannel || result instanceof DataflowBroadcast ) {
            inChannel = CH.getReadChannel(result)
        }
        else {
            inChannel = CH.value()
            inChannel.bind(result)
        }

        return inChannel
    }

    DataflowReadChannel getChannel() {
        return channel
    }

}
