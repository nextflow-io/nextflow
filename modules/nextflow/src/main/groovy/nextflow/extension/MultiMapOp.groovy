/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Channel
import nextflow.script.ChannelOut
import nextflow.script.TokenMultiMapDef
/**
 * Implements multi-map operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MultiMapOp {

    private DataflowReadChannel source
    private TokenMultiMapDef forkDef
    private Map<String, DataflowWriteChannel> targets = new LinkedHashMap<>(10)
    private ChannelOut output

    MultiMapOp(DataflowReadChannel source, Closure<TokenMultiMapDef> action) {
        this.source = source
        this.forkDef = action.call()

        for( String name : forkDef.names ) {
            targets.put( name, CH.createBy(source) )
        }
        this.output = new ChannelOut(targets)
    }

    ChannelOut getOutput() { this.output }

    protected void doNext(it) {
        final ret = (Map<String,?>)forkDef.closure.call(it)
        for( Map.Entry<String,?> entry : ret.entrySet() ) {
            targets[entry.key].bind(entry.value)
        }
    }

    protected void doComplete(nope) {
        for( DataflowWriteChannel ch : targets.values() ) {
            if( ch instanceof DataflowExpression ) {
                if( !ch.isBound()) ch.bind(Channel.STOP)
            }
            else {
                ch.bind(Channel.STOP)
            }
        }
    }

    MultiMapOp apply() {
        def events = new HashMap<String,Closure>(2)
        events.put('onNext', this.&doNext)
        events.put('onComplete', this.&doComplete)
        DataflowHelper.subscribeImpl(source, events)
        return this
    }
}
