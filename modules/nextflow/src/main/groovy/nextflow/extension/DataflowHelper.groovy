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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowEventListener
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import static java.util.Arrays.asList
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
        boolean onException(final DataflowProcessor processor, final Throwable e) {
            OperatorEx.log.error("@unknown", e)
            session.abort(e)
            return true;
        }
    }

    @PackageScope
    static DataflowEventAdapter stopErrorListener(DataflowReadChannel source, DataflowWriteChannel target) {

        new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( source instanceof DataflowExpression ) {
                    if( !(target instanceof DataflowExpression) )
                        processor.bindOutput( Channel.STOP )
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
    @PackageScope
    static DataflowProcessor newOperator( Map params, Closure code ) {

        // -- add a default error listener
        if( !params.containsKey('listeners') ) {
            // add the default error handler
            params.listeners = [ DEF_ERROR_LISTENER ]
        }

        final op = Dataflow.operator(params, code)
        NodeMarker.appendOperator(op)
        if( session && session.allOperators != null ) {
            session.allOperators.add(op)
        }

        return op
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
    @PackageScope
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
    @PackageScope
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
    @PackageScope
    static DataflowProcessor newOperator( DataflowReadChannel input, DataflowWriteChannel output, DataflowEventListener listener, Closure code ) {

        if( !listener )
            listener = DEF_ERROR_LISTENER

        def params = [:]
        params.inputs = [input]
        params.outputs = [output]
        params.listeners = [listener]

        final op = Dataflow.operator(params, code)
        NodeMarker.appendOperator(op)
        if( session && session.allOperators != null ) {
            session.allOperators << op
        }
        return op
    }

    /*
     * the list of valid subscription handlers
     */
    static private VALID_HANDLERS = [ 'onNext', 'onComplete', 'onError' ]

    /**
     * Verify that the map contains only valid names of subscribe handlers.
     * Throws an {@code IllegalArgumentException} when an invalid name is specified
     *
     * @param handlers The handlers map
     */
    @PackageScope
    static checkSubscribeHandlers( Map handlers ) {

        if( !handlers ) {
            throw new IllegalArgumentException("You must specify at least an event between: onNext, onComplete, onError")
        }

        handlers.keySet().each {
            if( !VALID_HANDLERS.contains(it) )  throw new IllegalArgumentException("Not a valid handler name: $it")
        }

    }

    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    static final DataflowProcessor subscribeImpl(final DataflowReadChannel source, final Map<String,Closure> events ) {
        checkSubscribeHandlers(events)

        def error = false
        def stopOnFirst = source instanceof DataflowExpression
        def listener = new DataflowEventAdapter() {

            @Override
            void afterStop(final DataflowProcessor processor) {
                if( !events.onComplete || error ) return
                try {
                    events.onComplete.call(processor)
                }
                catch( Exception e ) {
                    OperatorEx.log.error("@unknown", e)
                    session.abort(e)
                }
            }

            @Override
            boolean onException(final DataflowProcessor processor, final Throwable e) {
                error = true
                if( !events.onError ) {
                    log.error("@unknown", e)
                    session.abort(e)
                }
                else {
                    events.onError.call(e)
                }
                return true
            }
        }


        final Map<String, Object> parameters = new HashMap<String, Object>();
        parameters.put("inputs", [source])
        parameters.put("outputs", [])
        parameters.put('listeners', [listener])

        newOperator (parameters) {
            if( events.onNext ) {
                events.onNext.call(it)
            }
            if( stopOnFirst ) {
                ((DataflowProcessor) getDelegate()).terminate()
            }
        }
    }


    static DataflowProcessor chainImpl(final DataflowReadChannel source, final DataflowWriteChannel target, final Map params, final Closure closure) {

        final Map<String, Object> parameters = new HashMap<String, Object>(params)
        parameters.put("inputs", asList(source))
        parameters.put("outputs", asList(target))

        newOperator(parameters, new ChainWithClosure(closure))
    }

    /**
     * Implements the {@code #reduce} operator
     *
     * @param channel
     * @param seed
     * @param closure
     * @return
     */
    static DataflowProcessor reduceImpl(final DataflowReadChannel channel, final DataflowVariable result, def seed, final Closure closure) {

        // the *accumulator* value
        def accum = seed

        // intercepts operator events
        def listener = new DataflowEventAdapter() {
            /*
             * call the passed closure each time
             */
            void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                final item = messages.get(0)
                final value = accum == null ? item : closure.call(accum, item)

                if( value == Channel.VOID ) {
                    // do nothing
                }
                else if( value == Channel.STOP ) {
                    processor.terminate()
                }
                else {
                    accum = value
                }
            }

            /*
             * when terminates bind the result value
             */
            void afterStop(final DataflowProcessor processor) {
                result.bind(accum)
            }

            boolean onException(final DataflowProcessor processor, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        chainImpl(channel, CH.create(), [listeners: [listener]], {true})
    }

    @PackageScope
    @CompileStatic
    static KeyPair makeKey(List<Integer> pivot, entry) {
        final result = new KeyPair()

        if( !(entry instanceof List) ) {
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
                result.addValue(list[i])
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
