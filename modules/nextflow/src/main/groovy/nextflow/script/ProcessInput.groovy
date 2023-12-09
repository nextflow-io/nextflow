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
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.extension.CH
import nextflow.extension.ToListOp

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

    /**
     * Flag to support `each` input
     */
    private boolean iterator

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

        def value = obj instanceof Closure
            ? obj.call()
            : obj

        if( value == null )
            throw new IllegalArgumentException('A process input channel evaluates to null')

        if( iterator )
            value = getIteratorChannel(value)

        if( value instanceof DataflowReadChannel || value instanceof DataflowBroadcast )
            return CH.getReadChannel(value)

        final result = CH.value()
        result.bind(value)
        return result
    }

    private DataflowReadChannel getIteratorChannel(Object value) {
        def result
        if( value instanceof DataflowExpression ) {
            result = value
        }
        else if( CH.isChannel(value) ) {
            def read = CH.getReadChannel(value)
            result = new ToListOp(read).apply()
        }
        else {
            result = new DataflowVariable()
            result.bind(value)
        }

        return result.chainWith { it instanceof Collection || it == null ? it : [it] }
    }

    DataflowReadChannel getChannel() {
        return channel
    }

    void setIterator(boolean iterator) {
        this.iterator = iterator
    }

    boolean isIterator() {
        return iterator
    }

}
