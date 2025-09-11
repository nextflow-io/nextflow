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
 */

package nextflow.extension

import java.lang.reflect.InvocationTargetException

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowEventListener
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.Op
import nextflow.prov.OperatorRun
/**
 * This class provides helper methods to implement nextflow operators
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DataflowHelper {

    private static Session getSession() { Global.getSession() as Session }

    /**
     * Check if a {@code DataflowProcessor} is active
     *
     * @param operator A {@code DataflowProcessor} instance
     * @return {@code true} when the operator is still active or {@code false} otherwise eg. it received a poison pill
     */
    static boolean isProcessorActive( DataflowProcessor operator ) {
        def clazz = operator.class
        while( clazz != DataflowProcessor )
            clazz = operator.class.superclass
        def field = clazz.getDeclaredField('actor')
        field.setAccessible(true)
        def actor = field.get(operator)
        def method = actor.getClass().getMethod('isActive')
        return method.invoke(actor)
    }

    /*
     * The default operators listener when no other else is specified
     */
    static DataflowEventListener defaultErrorListener() {
        return new DataflowEventAdapter() {
            @Override
            boolean onException(final DataflowProcessor dp, final Throwable t) {
                final e = t instanceof InvocationTargetException ? t.cause : t
                OperatorImpl.log.error("@unknown", e)
                session?.abort(e)
                return true;
            }
        }
    }

    static DataflowEventListener stopErrorListener(DataflowReadChannel source, DataflowWriteChannel target) {
        new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor dp, final List<Object> messages) {
                if( source instanceof DataflowExpression ) {
                    if( target !instanceof DataflowExpression )
                        Op.bind(dp, target, Channel.STOP )
                    dp.terminate()
                }
            }

            @Override
            boolean onException(final DataflowProcessor dp, final Throwable e) {
                DataflowHelper.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }
    }

    static final DataflowProcessor subscribeImpl(final DataflowReadChannel source, final Map<String,Closure> events ) {
        new SubscribeOp()
            .withInput(source)
            .withEvents(events)
            .apply()
    }

    @PackageScope
    @CompileStatic
    static KeyPair makeKey(List<Integer> pivot, entry, OperatorRun run) {
        final result = new KeyPair()

        if( entry !instanceof List ) {
            if( pivot != [0] )
                throw new IllegalArgumentException("Not a valid `by` index: $pivot")
            result.keys = [entry]
            result.values = []
            return result
        }

        def list = (List)entry
        result.keys = new ArrayList(pivot.size())
        result.values = new ArrayList(list.size())

        for( int i=0; i<list.size(); i++ ) {
            if( i in pivot )
                result.addKey(list[i])
            else
                result.addValue(list[i], run)
        }

        return result
    }

    @PackageScope
    @CompileStatic
    static void addToList(List result, entry)  {
        if( entry instanceof List ) {
            result.addAll(entry)
        }
        else {
            result.add(entry)
        }
    }

}
