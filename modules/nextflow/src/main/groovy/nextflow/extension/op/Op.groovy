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

package nextflow.extension.op

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Global
import nextflow.Session
import nextflow.prov.OperatorRun
import nextflow.prov.Prov
import nextflow.prov.Tracker
/**
 * Operator helpers methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Op {

    static final public ConcurrentHashMap<DataflowProcessor, OpContext> context = new ConcurrentHashMap<>();

    static List<Object> unwrap(List messages) {
        final ArrayList<Object> result = new ArrayList<>();
        for( Object it : messages ) {
            result.add(it instanceof Tracker.Msg ? ((Tracker.Msg)it).value : it);
        }
        return result;
    }

    static Object unwrap(Object it) {
        return it instanceof Tracker.Msg ? it.value : it
    }

    static Tracker.Msg wrap(Object obj) {
        obj instanceof Tracker.Msg ? obj : Tracker.Msg.of(obj)
    }

    static void bind(DataflowProcessor operator, DataflowWriteChannel channel, List<Object> messages) {
        try {
            OperatorRun run=null
            for(Object msg : messages) {
                if( msg instanceof PoisonPill ) {
                    channel.bind(msg)
                    context.remove(operator)
                }
                else {
                    if( run==null ) {
                        final ctx = context.get(operator)
                        if( !ctx )
                            throw new IllegalStateException("Cannot find any context for operator=$operator")
                        run = ctx.getOperatorRun()
                    }
                    Prov.getTracker().bindOutput(run, channel, msg)
                }
            }
        }
        catch (Throwable t) {
            log.error("Unexpected resolving execution provenance: ${t.message}", t)
            (Global.session as Session).abort(t)
        }
    }


    static void bind(DataflowProcessor operator, DataflowWriteChannel channel, Object msg) {
        bind(operator, channel, List.of(msg))
    }

    static OpAbstractClosure instrument(Closure op, boolean accumulator=false) {
        return accumulator
            ? new OpGroupingClosure(op)
            : new OpRunningClosure(op)
    }

}
