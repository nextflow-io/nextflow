/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.DataflowEventListener
import nextflow.extension.op.Op
/**
 * Implements the chain operator
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ChainOp {

    private DataflowReadChannel source
    private DataflowWriteChannel target
    private List<DataflowEventListener> listeners = List.of()
    private Closure action

    static ChainOp create() {
        new ChainOp()
    }

    ChainOp withSource(DataflowReadChannel source) {
        assert source
        this.source = source
        return this
    }

    ChainOp withTarget(DataflowWriteChannel target) {
        assert target
        this.target = target
        return this
    }

    ChainOp withListener(DataflowEventListener listener) {
        assert listener != null
        this.listeners = List.of(listener)
        return this
    }

    ChainOp withListeners(List<DataflowEventListener> listeners) {
        assert listeners != null
        this.listeners = listeners
        return this
    }

    ChainOp withAction(Closure action) {
        this.action = action
        return this
    }

    DataflowWriteChannel apply() {
        assert source
        assert target
        assert action

        new Op()
            .withInput(source)
            .withOutput(target)
            .withListeners(listeners)
            .withCode(new ChainWithClosure(action))
            .apply()

        return target
    }
}
