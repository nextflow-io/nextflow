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
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
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
     * Create a dataflow object by the type of the specified source argument
     *
     * @param source
     * @return
     */
    @Deprecated
    static <V> DataflowChannel<V> newChannelBy(DataflowReadChannel<?> source) {

        switch( source ) {
            case DataflowExpression:
                return new DataflowVariable<V>()

            case DataflowQueue:
                return new DataflowQueue<V>()

            default:
                throw new IllegalArgumentException()
        }

    }


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
    @PackageScope
    static DEF_ERROR_LISTENER = new DataflowEventAdapter() {
        @Override
        boolean onException(final DataflowProcessor processor, final Throwable t) {
            final e = t instanceof InvocationTargetException ? t.cause : t
            OperatorImpl.log.error("@unknown", e)
            session?.abort(e)
            return true;
        }
    }

    @PackageScope
    static DataflowEventAdapter stopErrorListener(DataflowReadChannel source, DataflowWriteChannel target) {

        new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( source instanceof DataflowExpression ) {
                    if( target !instanceof DataflowExpression )
                        Op.bind(processor, target, Channel.STOP )
                    processor.terminate()
                }
            }

            @Override
            boolean onException(final DataflowProcessor processor, final Throwable e) {
                DataflowHelper.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

    }

    @PackageScope
    static Map createOpParams(inputs, outputs, listeners) {
        final params = new HashMap(3)
        params.inputs = inputs instanceof List ? inputs : [inputs]
        params.outputs = outputs instanceof List ? outputs : [outputs]
        params.listeners = listeners instanceof List ? listeners : [listeners]
        return params
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see nextflow.Session#allOperators
     *
     * @param params The map holding inputs, outputs channels and other parameters
     * @param code The closure to be executed by the operator
     */
    @Deprecated
    static DataflowProcessor newOperator( Map params, Closure code ) {

        // -- add a default error listener
        if( !params.listeners ) {
            // add the default error handler
            params.listeners = [ DEF_ERROR_LISTENER ]
        }

        return new Op()
            .withParams(params)
            .withCode(code)
            .apply()
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see nextflow.Session#allOperators
     *
     * @param inputs The list of the input {@code DataflowReadChannel}s
     * @param outputs The list of list output {@code DataflowWriteChannel}s
     * @param code The closure to be executed by the operator
     */
    @Deprecated
    static DataflowProcessor newOperator( List inputs, List outputs, Closure code ) {
        newOperator( inputs: inputs, outputs: outputs, code )
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see nextflow.Session#allOperators
     *
     * @param input An instance of {@code DataflowReadChannel} representing the input channel
     * @param output An instance of {@code DataflowWriteChannel} representing the output channel
     * @param code The closure to be executed by the operator
     */
    static DataflowProcessor newOperator( DataflowReadChannel input, DataflowWriteChannel output, Closure code ) {
        newOperator(input, output, DEF_ERROR_LISTENER, code )
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see nextflow.Session#allOperators
     *
     * @param input An instance of {@code DataflowReadChannel} representing the input channel
     * @param output An instance of {@code DataflowWriteChannel} representing the output channel
     * @param listener An instance of {@code DataflowEventListener} listening to operator's events
     * @param code The closure to be executed by the operator
     */
    static DataflowProcessor newOperator( DataflowReadChannel input, DataflowWriteChannel output, DataflowEventListener listener, Closure code ) {
        if( !listener )
            listener = DEF_ERROR_LISTENER

        final params = [:]
        params.inputs = [input]
        params.outputs = [output]
        params.listeners = [listener]

        return new Op()
            .withParams(params)
            .withCode(code)
            .apply()
    }

    static final DataflowProcessor subscribeImpl(final DataflowReadChannel source, final Map<String,Closure> events ) {
        new SubscribeOp()
            .withSource(source)
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

    @CompileStatic
    static Map<String,Closure> eventsMap(Closure onNext, Closure onComplete) {
        def result = new HashMap<String,Closure>(2)
        result.put('onNext', onNext)
        result.put('onComplete', onComplete)
        return result
    }
}
