/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.extension
import static java.util.Arrays.asList
import static nextflow.util.CheckHelper.checkParams

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowEventListener
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
/**
 * A set of operators inspired to RxJava extending the methods available on DataflowChannel
 * data structure
 *
 * See https://github.com/Netflix/RxJava/wiki/Observable
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class DataflowExtensions {

    private static Session getSession() { Global.getSession() as Session }

    /**
     * INTERNAL ONLY API
     * <p>
     * Add the {@code update} method to an {@code Agent} so that it call implicitly
     * the {@code Agent#updateValue} method
     *
     */
    static void update( Agent self, Closure message ) {
        assert message != null

        self.send {
            message.call(it)
            updateValue(it)
        }

    }

    /*
     * The default operators listener when no other else is specified
     */
    @PackageScope
    static DEF_ERROR_LISTENER = new DataflowEventAdapter() {
        @Override
        public boolean onException(final DataflowProcessor processor, final Throwable e) {
            DataflowExtensions.log.error("@unknown", e)
            session.abort(e)
            return true;
        }
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see Session#allProcessors
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
        if( session && session.allProcessors != null ) {
            session.allProcessors.add(op)
        }

        return op
    }

    /**
     * Creates a new {@code Dataflow.operator} adding the created instance to the current session list
     *
     * @see Session#allProcessors
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
     * @see Session#allProcessors
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
     * @see Session#allProcessors
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
        if( session && session.allProcessors != null ) {
            session.allProcessors << op
        }
        return op
    }


    /**
     * Create a dataflow object by the type of the specified source argument
     *
     * @param source
     * @return
     */
    @PackageScope
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
    static private final checkSubscribeHandlers( Map handlers ) {

        if( !handlers ) {
            throw new IllegalArgumentException("You must specify at least an event between: onNext, onComplete, onError")
        }

        handlers.keySet().each {
            if( !VALID_HANDLERS.contains(it) )  throw new IllegalArgumentException("Not a valid handler name: $it")
        }

    }

    /**
     * Subscribe *onNext* event
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> subscribe(final DataflowReadChannel<V> source, final Closure<V> closure) {
        subscribeImpl( source, [onNext: closure] )
        session.dag.addOperatorNode('subscribe', source, null)
        return source
    }

    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> subscribe(final DataflowReadChannel<V> source, final Map<String,Closure> events ) {
        subscribeImpl(source, events)
        session.dag.addOperatorNode('subscribe', source, null)
        return source
    }

    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowProcessor subscribeImpl(final DataflowReadChannel<V> source, final Map<String,Closure> events ) {
        checkSubscribeHandlers(events)

        def error = false
        def stopOnFirst = source instanceof DataflowExpression
        def listener = new DataflowEventAdapter() {

            @Override
            public void afterStop(final DataflowProcessor processor) {
                if( !events.onComplete || error ) return
                events.onComplete.call(processor)
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                error = true
                if( !events.onError ) {
                    DataflowExtensions.log.error("@unknown", e)
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

        newOperator(parameters) {
            if( events.onNext ) {
                events.onNext.call(it)
            }
            if( stopOnFirst ) {
                ((DataflowProcessor) getDelegate()).terminate()
            }
        }
    }

    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    public static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Closure<V> closure) {

        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>(closure))
        session.dag.addOperatorNode('chain', source, target)

        return target;
    }


    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    public static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Map<String, Object> params, final Closure<V> closure) {

        final DataflowReadChannel<V> target = newChannelBy(source)
        chainImpl(source, target, params, closure)
        session.dag.addOperatorNode('chain', source, target)

        return target;
    }

    private static <V> DataflowProcessor chainImpl(final DataflowReadChannel<?> source, final DataflowReadChannel<V> target, final Map params, final Closure<V> closure) {

        final Map<String, Object> parameters = new HashMap<String, Object>(params)
        parameters.put("inputs", asList(source))
        parameters.put("outputs", asList(target))

        newOperator(parameters, new ChainWithClosure<V>(closure))
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> map(final DataflowReadChannel<?> source, final Closure<V> closure) {
        assert source != null
        assert closure

        def target = new MapOp(source, closure).apply()

        session.dag.addOperatorNode('map', source, target)
        return target;
    }


    /**
     * Transform the items emitted by a channel by applying a function to each of them and then flattens the results of that function.
     *
     * @param source The source channel
     * @param closure The closure mapping the values emitted by the source channel
     * @return The channel emitting the mapped values
     */
    static public final <V> DataflowReadChannel<V> flatMap(final DataflowReadChannel<?> source, final Closure<V> closure=null) {
        assert source != null

        final target = new DataflowQueue()

        def listener = new DataflowEventAdapter() {
            @Override
            public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( source instanceof DataflowExpression ) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        newOperator(source, target, listener) {  item ->

            def result = closure != null ? closure.call(item) : item
            def proc = ((DataflowProcessor) getDelegate())

            switch( result ) {
                case Collection:
                    result.each { it -> proc.bindOutput(it) }
                    break

                case (Object[]):
                    result.each { it -> proc.bindOutput(it) }
                    break

                case Map:
                    result.each { it -> proc.bindOutput(it) }
                    break

                case Map.Entry:
                    proc.bindOutput( (result as Map.Entry).key )
                    proc.bindOutput( (result as Map.Entry).value )
                    break

                case Channel.VOID:
                    break

                default:
                    proc.bindOutput(result)
            }
        }

        session.dag.addOperatorNode('flatMap', source, target)
        return target
    }

    /**
     *
     * The reduce( ) operator applies a function of your choosing to the first item emitted by a source channel,
     * then feeds the result of that function along with the second item emitted by the source channel into the same
     * function, then feeds the result of that function along with the third item into the same function, and so on until
     * all items have been emitted by the source channel.
     *
     * Finally it emits the final result from the final call to your function as the sole output from the returned channel.
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> reduce(final DataflowReadChannel<?> source, final Closure<V> closure) {
        assert source instanceof DataflowQueue
        final target = new DataflowVariable()
        reduceImpl( source, target, null, closure )
        session.dag.addOperatorNode('reduce', source, target)
        return target
    }


    /**
     *
     * The reduce( ) operator applies a function of your choosing to the first item emitted by a source channel,
     * then feeds the result of that function along with the second item emitted by the source channel into the same
     * function, then feeds the result of that function along with the third item into the same function, and so on until
     * all items have been emitted by the source channel.
     *
     * Finally it emits the final result from the final call to your function as the sole output from the returned channel.
     *
     * @param source
     * @parama seed
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> reduce(final DataflowReadChannel<?> source, V seed, final Closure<V> closure) {
        assert !(source instanceof DataflowExpression)
        final target = new DataflowVariable()
        reduceImpl( source, target, seed, closure )
        session.dag.addOperatorNode('reduce', source, target)
        return target
    }

    /**
     * Implements the {@code #reduce} operator
     *
     * @param channel
     * @param seed
     * @param closure
     * @return
     */
    @PackageScope
    static <V> DataflowProcessor reduceImpl(final DataflowReadChannel<?> channel, final DataflowVariable result, def seed, final Closure<V> closure) {

        // the *accumulator* value
        def accum = seed

        // intercepts operator events
        def listener = new DataflowEventAdapter() {
            /*
             * call the passed closure each time
             */
            public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
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
            public void afterStop(final DataflowProcessor processor) {
                result.bind(accum)
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        chainImpl(channel, new DataflowQueue(), [listeners: [listener]], {true})
    }

    static public final DataflowReadChannel collectFile( final DataflowReadChannel source, final Closure closure = null ) {
        final result = new CollectFileOp(source, null, closure).apply()
        session.dag.addOperatorNode('collectFile', source, result)
        return result
    }

    static public final DataflowReadChannel collectFile( final DataflowReadChannel source, Map params, final Closure closure = null ) {
        def result = new CollectFileOp(source, params, closure).apply()
        session.dag.addOperatorNode('collectFile', source, result)
        return result
    }

    static public final DataflowReadChannel groupTuple( final DataflowReadChannel source, final Map params ) {
        def result = new GroupTupleOp(params, source).apply()
        session.dag.addOperatorNode('groupTuple', source, result)
        return result
    }

    /**
     * Iterates over the collection of items and returns each item that matches the given filter
     * by calling the {@code Object#isCase}method used by switch statements.
     *
     * This method can be used with different kinds of filters like regular expressions, classes, ranges etc. Example:
     *
     * def list = ['a', 'b', 'aa', 'bc', 3, 4.5]
     * assert list.filter( ~/a+/ )  == ['a', 'aa']
     * assert list.filter( ~/../ )  == ['aa', 'bc']
     * assert list.filter( Number ) == [ 3, 4.5 ]
     * assert list.filter{ it.toString().size() == 1 } == [ 'a', 'b', 3 ]
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final <V> DataflowReadChannel<V> filter(final DataflowReadChannel<V> source, final Object criteria) {
        def discriminator = new BooleanReturningMethodInvoker("isCase");
        def target = newChannelBy(source)
        newOperator(source, target, {
            def result = discriminator.invoke(criteria, (Object)it)
            if( result ) target.bind(it)
        })
        session.dag.addOperatorNode('filter', source, target)
        return target
    }

    static public <T> DataflowReadChannel<T> filter(DataflowReadChannel<T> source, final Closure<Boolean> closure) {
        def target = newChannelBy(source)
        newOperator(source, target, {
            def result = DefaultTypeTransformation.castToBoolean(closure.call(it))
            if( result ) target.bind(it)
        })
        session.dag.addOperatorNode('filter', source, target)
        return target
    }

    static public <T> DataflowReadChannel<T> until(DataflowReadChannel<T> source, final Closure<Boolean> closure) {
        def target = newChannelBy(source)
        newOperator(source, target, {
            final result = DefaultTypeTransformation.castToBoolean(closure.call(it))
            final proc = ((DataflowProcessor) getDelegate())

            if( result ) {
                proc.bindOutput(Channel.STOP)
                proc.terminate()
            }
            else {
                proc.bindOutput(it)
            }
        })
        session.dag.addOperatorNode('until', source, target)
        return target
    }


    /**
     * Modifies this collection to remove all duplicated items, using the default comparator.
     *
     * assert [1,3] == [1,3,3].unique()
     *
     * @param source
     * @return
     */
    static public final <V> DataflowReadChannel<V> unique(final DataflowReadChannel<V> source) {
        unique(source) { it }
    }

    /**
     * A convenience method for making a collection unique using a Closure to determine duplicate (equal) items. If the closure takes a single parameter, the argument passed will be each element, and the closure should return a value used for comparison (either using Comparable#compareTo or Object#equals). If the closure takes two parameters, two items from the collection will be passed as arguments, and the closure should return an int value (with 0 indicating the items are not unique).
     * assert [1,4] == [1,3,4,5].unique { it % 2 }
     * assert [2,3,4] == [2,3,3,4].unique { a, b -> a <=> b }
     *
     * @param source
     * @param comparator
     * @return
     */
    static public final <V> DataflowReadChannel<V> unique(final DataflowReadChannel<V> source, Closure comparator ) {

        def history = [:]
        def target = newChannelBy(source)

        // when the operator stop clear the history map
        def events = new DataflowEventAdapter() {
            public void afterStop(final DataflowProcessor processor) {
                history.clear()
                history = null
            }
        }

        def filter = {
            def key = comparator.call(it)
            if( history.containsKey(key) ) {
                return Channel.VOID
            }
            else {
                history.put(key,true)
                return it
            }
        }  as Closure<V>

        // filter removing all duplicates
        chainImpl(source, target, [listeners: [events]], filter )

        session.dag.addOperatorNode('unique', source, target)
        return target
    }

    /**
     * Suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     *
     * @return
     */
    static public final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> source ) {
        distinct(source) {it}
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    static public final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> source, Closure<?> comparator ) {

        def previous = null
        final DataflowReadChannel<V> target = newChannelBy(source)
        Closure<V> filter = { it ->

            def key = comparator.call(it)
            if( key == previous ) {
                return Channel.VOID
            }
            previous = key
            return it
        }

        chainImpl(source, target, [:], filter)

        session.dag.addOperatorNode('distinct', source, target)
        return target
    }

    /**
     *
     * Emit only the first item emitted by a channel, or the first item that meets some condition
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#first
     *
     * @param source
     * @return
     */
    static public final <V> DataflowReadChannel<V> first( DataflowReadChannel<V> source ) {

        def target = new DataflowVariable<V>()
        source.whenBound { target.bind(it) }
        session.dag.addOperatorNode('first', source, target)
        return target
    }

    /**
     *
     * Emit only the first item emitted by a channel, or the first item that meets some condition
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#first
     *
     * @param source
     * @return
     */
    static public final <V> DataflowReadChannel<V> first( final DataflowReadChannel<V> source, Object criteria ) {
        assert !(source instanceof DataflowExpression)

        def target = new DataflowVariable()
        def discriminator = new BooleanReturningMethodInvoker("isCase");

        newOperator([source],[]) {
            if( discriminator.invoke(criteria, it) ) {
                target.bind(it)
                ((DataflowProcessor) getDelegate()).terminate()
            }
        }

        session.dag.addOperatorNode('first', source, target)
        return target
    }

    /**
     *
     * emit only the first n items emitted by an Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#take
     *
     * @param source
     * @param n The number of items to be taken. The value {@code -1} has a special semantic for all
     * @return The resulting channel emitting the taken values
     */
    static public final <V> DataflowReadChannel<V> take( final DataflowReadChannel<V> source, int n ) {
        assert !(source instanceof DataflowExpression)

        def count = 0
        final target = new DataflowQueue<V>()

        if( n==0 ) {
            target.bind(Channel.STOP)
            return target
        }

        final listener = new DataflowEventAdapter() {
            @Override
            public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( ++count >= n ) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        newOperator(
                inputs: [source],
                outputs: [target],
                listeners: (n > 0 ? [listener] : []),
                new ChainWithClosure(new CopyChannelsClosure()))

        session.dag.addOperatorNode('take', source, target)
        return target
    }

    /**
     * The last operator creates a channel that only returns the last item emitted by the source channel
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the `last` item in the channel
     */
    static public final <V> DataflowReadChannel<V> last( final DataflowReadChannel<V> source ) {

        def target = new DataflowVariable()
        def V last = null
        subscribeImpl( source, [onNext: { last = it }, onComplete: {  target.bind(last) }] )
        session.dag.addOperatorNode('last', source, target)
        return target
    }


    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static public final <V> DataflowReadChannel<V> toList(final DataflowReadChannel<V> source) {
        final target = ToListOp.apply(source)
        session.dag.addOperatorNode('toList', source, target)
        return target
    }

    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List} sorting its content
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static public final <V> DataflowReadChannel<V> toSortedList(final DataflowReadChannel<V> source, Closure closure = null) {
        final reduced = new DataflowVariable()
        reduceImpl(source, reduced, []) { list, item -> list << item }
        final target = reduced.then { List list ->
            closure ? list.sort(closure) : list.sort()
        }
        session.dag.addOperatorNode('toSortedList', source, target)
        return target as DataflowVariable
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param source
     * @param value
     * @return
     */
    static public final DataflowReadChannel<Number> count(final DataflowReadChannel<?> source ) {
        final target = count0(source, null)
        session.dag.addOperatorNode('count', source, target)
        return target
    }

    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param source
     * @param criteria
     * @return
     */
    static public final DataflowReadChannel<Number> count(final DataflowReadChannel<?> source, final Object criteria ) {
        final target = count0(source, criteria)
        session.dag.addOperatorNode('count', source, target)
        return target
    }

    private static DataflowVariable count0(DataflowReadChannel<?> source, Object criteria) {

        final target = new DataflowVariable()
        final discriminator = criteria != null ? new BooleanReturningMethodInvoker("isCase") : null

        if( source instanceof DataflowExpression) {
            source.whenBound { item ->
                discriminator == null || discriminator.invoke(criteria, item) ? target.bind(1) : target.bind(0)
            }
        }
        else {
            reduceImpl(source, target, 0) { current, item ->
                discriminator == null || discriminator.invoke(criteria, item) ? current+1 : current
            }
        }

        return target
    }


    /**
     * Groups the items emitted by the source channel into groups determined by the supplied mapping closure and counts the frequency of the created groups
     * @param source The source channel
     * @return A {@code DataflowVariable} returning the a {@code Map} containing the counting values for each key
     */
    static public final DataflowReadChannel<Map> countBy(final DataflowReadChannel<?> source ) {
        countBy(source, { it })
    }

    /**
     * Sorts all collection members into groups determined by the supplied mapping closure and counts the group size
     *
     * @param source
     * @param criteria
     * @return
     */
    static public final DataflowReadChannel<Map> countBy(final DataflowReadChannel<?> source, final Closure criteria ) {

        final target = new DataflowVariable()

        reduceImpl(source, target, [:]) { Map map, item ->
                def key = criteria.call(item)
                def value = map.containsKey(key) ? map.get(key)+1 : 1
                map.put(key, value)
                return map
        }

        session.dag.addOperatorNode('countBy', source, target)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} returning the minimum value
     */
    static public final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> source) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { min, val -> val<min ? val : min }
        session.dag.addOperatorNode('min', source, target)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @param comparator If the closure has two parameters it is used like a traditional Comparator. I.e. it should compare
     *      its two parameters for order, returning a negative integer, zero, or a positive integer when the first parameter
     *      is less than, equal to, or greater than the second respectively. Otherwise, the Closure is assumed to take a single
     *      parameter and return a Comparable (typically an Integer) which is then used for further comparison.
     * @return  A {@code DataflowVariable} returning the minimum value
     */
    static public final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> source, Closure<V> comparator) {

        def action
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure<V>){ min, item -> comparator.call(item) < comparator.call(min) ? item : min  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure<V>){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        final target = new DataflowVariable()
        reduceImpl(source, target, null, action)
        session.dag.addOperatorNode('min', source, target)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @param comparator The a {@code Comparator} object
     * @return A {@code DataflowVariable} returning the minimum value
     */
    static public final <V> DataflowReadChannel<V> min(final DataflowQueue<V> source, Comparator comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)<0 ? a : b }
        session.dag.addOperatorNode('min', source, target)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    static public final <V> DataflowReadChannel<V> max(final DataflowQueue source) {
        final target = new DataflowVariable()
        reduceImpl(source,target, null) { max, val -> val>max ? val : max }
        session.dag.addOperatorNode('max', source, target)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @param comparator If the closure has two parameters it is used like a traditional Comparator. I.e. it should compare
     *  its two parameters for order, returning a negative integer, zero, or a positive integer when the first parameter is
     *  less than, equal to, or greater than the second respectively. Otherwise, the Closure is assumed to take a single
     *  parameter and return a Comparable (typically an Integer) which is then used for further comparison
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    static public final <V> DataflowReadChannel<V> max(final DataflowQueue<V> source, Closure comparator) {

        def action
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure<V>){ max, item -> comparator.call(item) > comparator.call(max) ? item : max  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure<V>){ a, b ->  comparator.call(a,b)>0 ? a : b  }
        }
        else {
            throw new IllegalArgumentException("Comparator closure can accept at most 2 arguments")
        }

        final target = new DataflowVariable()
        reduceImpl(source, target, null, action)
        session.dag.addOperatorNode('max', source, target)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @param comparator A {@code Comparator} object
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    static public final <V> DataflowVariable<V> max(final DataflowQueue<V> source, Comparator<V> comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)>0 ? a : b }
        session.dag.addOperatorNode('max', source, target)
        return target
    }

    /**
     * The sum operators crates a channel that emits the sum of all values emitted by the source channel to which is applied
     *
     * @param source  The source channel providing the values to sum
     * @param closure  A closure that given an entry returns the value to sum
     * @return A {@code DataflowVariable} emitting the final sum value
     */
    static public final DataflowReadChannel sum(final DataflowQueue source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'sum', action: closure)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
        session.dag.addOperatorNode('sum', source, target)
        return target
    }


    static public final DataflowReadChannel mean(final DataflowQueue source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'mean', action: closure, mean: true)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
        session.dag.addOperatorNode('mean', source, target)
        return target
    }

    private static class Aggregate {

        def accum
        long count = 0
        boolean mean
        Closure action
        String name

        def process(it) {
            if( it == null || it == Channel.VOID )
                return

            count++

            def item = action != null ? action.call(it) : it
            if( accum == null )
                accum = item

            else if( accum instanceof Number )
                accum += item

            else if( accum instanceof List && item instanceof List)
                for( int i=0; i<accum.size() && i<item.size(); i++ )
                    accum[i] += item.get(i)

            else
                throw new IllegalArgumentException("Invalid `$name` item: $item [${item.class.simpleName}]")
        }

        def getResult() {
            if( !mean || count == 0 )
                return accum

            if( accum instanceof List )
                return accum.collect { it / count }
            else
                return accum / count
        }
    }

    /**
     * Sorts all collection members into groups determined by the supplied mapping closure
     *
     * @param source
     * @param mapper
     * @return
     */
    static public final DataflowReadChannel<Map> groupBy(final DataflowReadChannel source, final params = null ) {

        int index = 0
        Closure mapper = DEFAULT_MAPPING_CLOSURE

        if( params instanceof Closure )
            mapper = params

        else if( params instanceof Number ) {
            index = params as int
        }
        else if( params != null ) {
            throw new IllegalArgumentException("Not a valid `group` argument: $params")
        }

        final target = new DataflowVariable()
        final int len = mapper.getMaximumNumberOfParameters()
        reduceImpl(source, target, [:]) { map, item ->
            def key = len == 2 ? mapper.call(item,index) : mapper.call(item)
            def list = map.get(key)
            list = list ? list << item : [item]
            map.put(key, list)
            return map
        }

        session.dag.addOperatorNode('groupBy', source, target)
        return target
    }


    /**
     * Given a an associative array mapping a key with the destination channel, the operator route forwards the items emitted
     * by the source channel to the target channel matching the key in the routing map
     *
     * @param source The source channel emitting the value to route
     * @param targets The routing map i.e. a {@code Map} associating each key to the target channel
     * @param mapper A optional mapping function that given an entry return its key
     */
    static public final void route( final DataflowReadChannel source, Map<?,DataflowWriteChannel> targets, Closure mapper = DEFAULT_MAPPING_CLOSURE ) {

        source.subscribeImpl (
                [
                        onNext: { value ->
                            def key = mapper ? mapper.call(value) : value
                            def channel = targets.get(key)
                            // emit the value itself
                            if( channel ) {
                                channel << value
                            }

                        },

                        onComplete: {
                            targets.values().each { it << Channel.STOP }
                        }

                ]
        )

        session.dag.addOperatorNode('route', source, targets.values())
    }

    static public final DataflowReadChannel route( final DataflowReadChannel source, final Closure mapper = DEFAULT_MAPPING_CLOSURE ) {
        assert !(source instanceof DataflowExpression)

        def allChannels = new ConcurrentHashMap()
        DataflowQueue target = new DataflowQueue()

        source.subscribeImpl (
                [
                    onNext: { value ->
                        def key = mapper ? mapper.call(value) : value
                        def channel = allChannels.get(key)
                        if( channel == null ) {
                            channel = new DataflowQueue()
                            allChannels[key] = channel
                            // emit the key - channel pair
                            target << [ key, channel ]
                        }
                        // emit the value itself
                        channel << value
                    },

                    onComplete: {
                        allChannels.values().each { it << Channel.STOP }
                        target << Channel.STOP
                    }

                ]
        )

        session.dag.addOperatorNode('route', source, target)
        return target
    }


    static public final DataflowReadChannel spread( final DataflowReadChannel source, Object other ) {

        final target = new DataflowQueue()

        def inputs
        switch(other) {
            case DataflowQueue: inputs = ToListOp.apply((DataflowQueue)other); break
            case DataflowExpression: inputs = other; break
            case Collection: inputs = Channel.value(other); break
            case (Object[]): inputs = Channel.value(other as List); break
            default: throw new IllegalArgumentException("Not a valid argument for 'spread' operator [${other?.class?.simpleName}]: ${other} -- Use a Collection or a channel instead. ")
        }

        final stopOnFirst = source instanceof DataflowExpression
        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst ) return
                processor.terminate()
                target.bind(Channel.STOP)
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        final params = [:]
        params.inputs = [source, inputs]
        params.outputs = [target]
        params.listeners = [listener]

        newOperator(params) { a, b ->
            def proc = ((DataflowProcessor) getDelegate())
            [ [a], (b as List) ]
                    .combinations()
                    .each{ Collection it -> proc.bindOutput(it.flatten())  }
        }

        def sources = new ArrayList(2)
        sources.add ( source )
        sources.add ( other instanceof DataflowChannel ? other : inputs  )
        session.dag.addOperatorNode('spread', sources, target)
        return target
    }


    static public final DataflowReadChannel flatten( final DataflowReadChannel source )  {

        final listeners = []
        final target = new DataflowQueue()

        if( source instanceof DataflowExpression ) {
            listeners << new DataflowEventAdapter() {
                @Override
                public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }

                public boolean onException(final DataflowProcessor processor, final Throwable e) {
                    DataflowExtensions.log.error("@unknown", e)
                    session.abort(e)
                    return true;
                }
            }
        }


        newOperator(inputs: [source], outputs: [target], listeners: listeners) {  item ->

            def proc = ((DataflowProcessor) getDelegate())
            switch( item ) {
                case Collection:
                    item.flatten().each { value -> proc.bindOutput(value) }
                    break

                case (Object[]):
                    item.flatten().each { value -> proc.bindOutput(value) }
                    break

                case Channel.VOID:
                    break

                default:
                    proc.bindOutput(item)
            }
        }

        session.dag.addOperatorNode('flatten', source, target)
        return target
    }

    /**
     * The ``buffer( )`` operator gathers the items emitted by the source channel into bundles and
     * and emits these bundles as its own emissions.
     *
     * @param source The dataflow channel from where the values are gathered
     * @param closingCriteria A condition that has to be verified to close
     * @return A newly created dataflow queue which emitted the gathered values as bundles
     */
    static public final <V> DataflowReadChannel<V> buffer( final DataflowReadChannel<V> source, Object closingCriteria ) {

        def target = new DataflowQueue()
        def closure = new BooleanReturningMethodInvoker("isCase");
        bufferImpl(source, target, null, { Object it -> closure.invoke(closingCriteria, it) }, false)

        session.dag.addOperatorNode('buffer', source, target)
        return target
    }

    static public final <V> DataflowReadChannel<V> buffer( final DataflowReadChannel<V> source, Object startingCriteria, Object closingCriteria ) {
        assert startingCriteria != null
        assert closingCriteria != null

        def c1 = new BooleanReturningMethodInvoker("isCase");
        def c2 = new BooleanReturningMethodInvoker("isCase");

        def target = new DataflowQueue()
        bufferImpl(source, target, {Object it -> c1.invoke(startingCriteria, it)}, {Object it -> c2.invoke(closingCriteria, it)}, false)

        session.dag.addOperatorNode('buffer', source, target)
        return target
    }

    static public final <V> DataflowReadChannel<V> buffer( DataflowReadChannel<V> source, Map<String,?> params ) {
        checkParams( 'buffer', params, 'size','skip','remainder' )

        int skip = (int)(params?.skip ?: 0)
        int size = params.size as int
        boolean remainder = params?.remainder ?: false
        if( !size ) {
            throw new IllegalArgumentException()
        }

        final target = new DataflowQueue()
        bufferWithSizeConstraint( source, target, size, skip, remainder )

        session.dag.addOperatorNode('buffer', source, target)
        return target
    }

    static private <V> DataflowProcessor bufferWithSizeConstraint( final DataflowReadChannel<V> source, DataflowQueue target, int size, int skip, boolean reminder ) {
        assert size>0

        int skipCount = 0
        int itemCount = 0

        def closeRule = {
            itemCount +=1
            if( itemCount-skip == size ) {
                itemCount = 0;
                return true
            }
            return false
        }


        def startRule = {
            skipCount +=1
            if( skipCount > skip ) {
                skipCount = 0
                return true
            }
            return false
        }

        bufferImpl(source, target, skip>0 ? startRule : null, closeRule, reminder )
    }


    static private <V> DataflowProcessor bufferImpl( DataflowReadChannel<V> source, DataflowQueue target, Closure startingCriteria, Closure closeCriteria, boolean remainder ) {
        assert closeCriteria

        // the list holding temporary collected elements
        def buffer = []

        // -- intercepts the PoisonPill and sent out the items remaining in the buffer when the 'remainder' flag is true
        def listener = new DataflowEventAdapter() {

            public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( message instanceof PoisonPill && remainder && buffer.size() )
                    target.bind(buffer)

                return message;
            }

            @Override
            boolean onException(DataflowProcessor processor, Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

        // -- open frame flag
        boolean isOpen = startingCriteria == null

        // -- the operator collecting the elements
        newOperator( source, target, listener ) {
            if( isOpen ) {
                buffer << it
            }
            else if( startingCriteria.call(it) ) {
                isOpen = true
                buffer << it
            }

            if( closeCriteria.call(it) ) {
                ((DataflowProcessor) getDelegate()).bindOutput(buffer);
                buffer = []
                // when a *startingCriteria* is defined, close the open frame flag
                isOpen = (startingCriteria == null)
            }
        }
    }

    static public final <V> DataflowReadChannel<V> collate( DataflowReadChannel<V> source, int size, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        final target = new DataflowQueue()
        bufferWithSizeConstraint( source, target, size, 0, keepRemainder )

        session.dag.addOperatorNode('collate', source, target)
        return target
    }

    static public final <V> DataflowReadChannel<V> collate( DataflowReadChannel<V> source, int size, int step, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        if( step <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'step' for operator 'collate' -- it must be greater than zero: $step")
        }

        // the result queue
        final target = new DataflowQueue();

        // the list holding temporary collected elements
        List<List<?>> allBuffers = []

        // -- intercepts the PoisonPill and sent out the items remaining in the buffer when the 'remainder' flag is true
        def listener = new DataflowEventAdapter() {

            public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( message instanceof PoisonPill && keepRemainder && allBuffers.size() ) {
                    allBuffers.each {
                        target.bind( it )
                    }
                }

                return message;
            }

            @Override
            boolean onException(DataflowProcessor processor, Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }


        int index = 0

        // -- the operator collecting the elements
        newOperator( inputs: [source], outputs: [target], listeners: [listener]) {

            if( index++ % step == 0 ) {
                allBuffers.add( [] )
            }

            allBuffers.each { List list -> list.add(it) }

            def buf = allBuffers.head()
            if( buf.size() == size )  {
                ((DataflowProcessor) getDelegate()).bindOutput(buf)
                allBuffers = allBuffers.tail()
            }

        }

        session.dag.addOperatorNode('collate', source, target)
        return target
    }


    /**
     * Similar to https://github.com/Netflix/RxJava/wiki/Combining-Observables#merge
     *
     * @param source
     * @param others
     * @return
     */
    static final DataflowReadChannel mix( DataflowReadChannel source, DataflowReadChannel... others ) {
        assert others.size()>0

        def target = new DataflowQueue()
        def count = new AtomicInteger( others.size()+1 )
        def handlers = [
                onNext: { target << it },
                onComplete: { if(count.decrementAndGet()==0) { target << Channel.STOP } }
        ]

        subscribeImpl(source, handlers)
        others.each{ it.subscribeImpl(handlers) }

        def allSources = [source]
        allSources.addAll(others)

        session.dag.addOperatorNode('mix', allSources, target)
        return target
    }


    /**
     * Phase channels
     *
     * @param source
     * @param other
     * @param mapper
     * @return
     */
    static final DataflowReadChannel phase( DataflowReadChannel source, Map opts, DataflowReadChannel other, Closure mapper = null ) {

        def target = new PhaseOp(source,other)
                        .setMapper(mapper)
                        .setOpts(opts)
                        .apply()

        session.dag.addOperatorNode('phase', [source, other], target)
        return target
    }

    static final DataflowReadChannel phase( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new PhaseOp(source,other)
                        .setMapper(mapper)
                        .apply()

        session.dag.addOperatorNode('phase', [source, other], target)
        return target
    }

    /**
     * Implements the default mapping strategy, having the following strategy:
     * <pre>
     *     Map -> first entry key
     *     Map.Entry -> the entry key
     *     Collection -> first item
     *     Array -> first item
     *     Object -> the object itself
     * </pre>
     * @param obj
     * @return
     */

    public static Closure DEFAULT_MAPPING_CLOSURE = { obj, int index=0 ->

        switch( obj ) {

            case List:
                def values = (List)obj
                return values.size() ? values.get(index) : null

            case (Object[]):
                def values = (Object[])obj
                return values.size() ? values[index] : null

            case Map:
                obj = ((Map)obj).values()
                // note: fallback into the following case

            case Collection:
                def itr = ((Collection)obj) .iterator()
                def count=0
                while( itr.hasNext() ) {
                    def value = itr.next()
                    if( count++ == index ) return value
                }
                return null

            case Map.Entry:
                def entry = (Map.Entry) obj
                return (index == 0 ? entry.key :
                        index == 1 ? entry.value : null)

            default:
                return index==0 ? obj : null
        }

    }


    public static <T> DataflowReadChannel cross( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new CrossOp(source, other)
                    .setMapper(mapper)
                    .apply()

        session.dag.addOperatorNode('cross', [source,other], target)
        return target
    }


    /**
     * Creates a channel that emits the items in same order as they are emitted by two or more channel
     *
     * @param source
     * @param others
     * @return
     */
    static final DataflowWriteChannel concat( DataflowReadChannel source, DataflowReadChannel... others ) {

        def target = new ConcatOp(source, others).apply()

        def allSources = [source]
        if(others) allSources.addAll(others)

        session.dag.addOperatorNode('concat', allSources, target)
        return target
    }


    /**
     * When the items emitted by the source channel are tuples of values, the operator separate allows you to specify a
     * list of channels as parameters, so that the value i-th in a tuple will be assigned to the target channel
     * with the corresponding position index.
     *
     * @param source The source channel
     * @param outputs An open array of target channels
     */
    static void separate( DataflowReadChannel source, final DataflowWriteChannel... outputs ) {
        new SeparateOp(source, outputs as List<DataflowWriteChannel>).apply()
        session.dag.addOperatorNode('separate', source, outputs)
    }

    static void separate(final DataflowReadChannel source, final List<DataflowWriteChannel<?>> outputs) {
        new SeparateOp(source, outputs).apply()
        session.dag.addOperatorNode('separate', source, outputs)
    }

    static void separate(final DataflowReadChannel source, final List<DataflowWriteChannel<?>> outputs, final Closure<List<Object>> code) {
        new SeparateOp(source, outputs, code).apply()
        session.dag.addOperatorNode('separate', source, outputs)
    }

    static public final List<DataflowReadChannel> separate( final DataflowReadChannel source, int n ) {
        def outputs = new SeparateOp(source, n).apply()
        session.dag.addOperatorNode('separate', source, outputs)
        return outputs
    }

    static public final List<DataflowReadChannel> separate( final DataflowReadChannel source, int n, Closure mapper  ) {
        def outputs = new SeparateOp(source, n, mapper).apply()
        session.dag.addOperatorNode('separate', source, outputs)
        return outputs
    }


    static final void into( DataflowReadChannel source, final DataflowWriteChannel... targets ) {
        new IntoOp(source, targets as List<DataflowWriteChannel>).apply()
        session.dag.addOperatorNode('into', source, targets)
    }

    static public final List<DataflowReadChannel> into( final DataflowReadChannel source, int n ) {
        def outputs = new IntoOp(source,n).apply().getOutputs()
        session.dag.addOperatorNode('into', source, outputs)
        return outputs
    }

    /**
     * Forward source dataflow channel *into* one or more dataflow channels. For example:
     * <pre>
     *     Channel.from( ... )
     *            .map { ... }
     *            .into { foo; bar }
     * </pre>
     *
     * It creates two new dataflow variables named {@code foo} and {@code bar} and copied the map
     * result into them.
     *
     * @param source The source dataflow channel which items are copied into newly created dataflow variables.
     * @param holder A closure that defines one or more variable names into which source items are copied.
     */
    static void into( DataflowReadChannel source, Closure holder ) {
        def outputs = new IntoOp(source,holder).apply().getOutputs()
        session.dag.addOperatorNode('into', source, outputs)
    }



    /**
     * Implements a tap that create implicitly a new dataflow variable in the global script context.
     * For example:
     *
     * <pre>
     *     Channel.from(...)
     *            .tap { newChannelName }
     *            .map { ... }
     *  </pre>
     *
     * @param source The source dataflow variable
     * @param holder The closure defining the new variable name
     * @return The tap resulting dataflow channel
     */
    static DataflowReadChannel tap( final DataflowReadChannel source, final Closure holder ) {
        def tap = new TapOp(source, holder).apply()
        session.dag.addOperatorNode('tap', source, tap.outputs)
        return (DataflowReadChannel)tap.result
    }

    static DataflowReadChannel tap( final DataflowReadChannel source, final DataflowWriteChannel target ) {
        def tap = new TapOp(source, target).apply()
        session.dag.addOperatorNode('tap', source, tap.outputs)
        return (DataflowReadChannel)tap.result
    }

    /**
     * Assign the {@code source} channel to a global variable with the name specified by the closure.
     * For example:
     * <pre>
     *     Channel.from( ... )
     *            .map { ... }
     *            .set { newChannelName }
     * </pre>
     *
     * @param DataflowReadChannel
     * @param holder A closure that must define a single variable expression
     */
    static void set( DataflowReadChannel source, Closure holder ) {
        final name = CaptureProperties.capture(holder)
        if( !name )
            throw new IllegalArgumentException("Missing name to which set the channel variable")

        if( name.size()>1 )
            throw new IllegalArgumentException("Operation `set` does not allow more than one target name")

        final binding = Global.session.binding
        binding.setVariable(name[0], source)
    }

    /**
     * Empty the specified value only if the source channel to which is applied is empty i.e. do not emit
     * any value.
     *
     * @param source The channel to which the operator is applied
     * @param value The value to emit when the source channel is empty. If a closure is used the the value returned by its invocation is used.
     * @return The resulting channel emitting the source items or the default value when the channel is empty
     */
    static public DataflowReadChannel ifEmpty( DataflowReadChannel source, value ) {

        boolean empty = true
        def result = newChannelBy(source)
        subscribeImpl(source, [
                onNext: { result.bind(it); empty=false },
                onComplete: {
                    if(empty) {
                        if( value instanceof Closure )
                            result.bind(value.call())
                        else
                            result.bind(value)
                    }
                    result.bind(Channel.STOP)
                }]
        )

        session.dag.addOperatorNode('ifEmpty', source, result)
        return result
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    static public void print(final DataflowReadChannel<?> source, Closure closure = null) {
        subscribeImpl(source, [onNext: { System.out.print( closure ? closure.call(it) : it ) }])
        session.dag.addOperatorNode('print', source, null)
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    static public void println(final DataflowReadChannel<?> source, Closure closure = null) {
        subscribeImpl(source, [onNext: { System.out.println( closure ? closure.call(it) : it ) }])
        session.dag.addOperatorNode('println', source, null)
    }


    static private final PARAMS_VIEW = [newLine: Boolean]

    /**
     * Print out the channel content retuning a new channel emitting the identical content as the original one
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> view(final DataflowReadChannel<?> source, Map opts, Closure closure = null) {
        assert source != null
        checkParams('view', opts, PARAMS_VIEW)
        final newLine = opts.newLine != false

        final target = newChannelBy(source);
        final isChannel = !(target instanceof DataflowExpression)

        final printHandle = newLine ? System.out.&println : System.out.&print

        final apply = [

                onNext:
                        {
                            printHandle ( closure != null ? closure.call(it) : it )
                            target.bind(it)
                        },

                onComplete: {
                    if(isChannel) target.bind(Channel.STOP)
                }
        ]

        subscribeImpl(source,apply)

        session.dag.addOperatorNode('view', source, target)
        return target;

    }

    static public final <V> DataflowReadChannel<V> view(final DataflowReadChannel<?> source, Closure closure = null) {
        view(source, [:], closure)
    }

    /**
     * Creates a channel emitting the entries in the collection to which is applied
     * @param values
     * @return
     */
    static DataflowQueue channel(Collection values) {
        def target = new DataflowQueue()
        def itr = values.iterator()
        while( itr.hasNext() ) target.bind(itr.next())
        target.bind(Channel.STOP)
        session.dag.addSourceNode('channel',target)
        return target
    }

    @Deprecated
    static public DataflowBroadcast broadcast( DataflowReadChannel source ) {
        log.warn("Operator `broadcast` is deprecated -- It will be removed in a future release")
        def result = new DataflowBroadcast()
        source.into(result)
        return result
    }

    /**
     * Close a dataflow queue channel binding a {@link Channel#STOP} item
     *
     * @param source The source dataflow channel to be closed.
     */
    static close( DataflowReadChannel source ) {
        if( source instanceof DataflowQueue ) {
            source.bind(Channel.STOP)
        }
        else {
            log.warn "Operation `close` can only be applied to a queue channel"
        }
        return source
    }

    static <T> void choice(final DataflowReadChannel source, final List<DataflowWriteChannel<T>> outputs, final Closure<Integer> code) {
        new ChoiceOp(source,outputs,code).apply()
        session.dag.addOperatorNode('choice', source, outputs)
    }

    static <V> DataflowReadChannel<V> merge(final DataflowReadChannel<V> source, final DataflowReadChannel<V> other, final Closure<V> closure) {
        final result = newChannelBy(source)
        final inputs = [source, other]
        newOperator(inputs, [result], new ChainWithClosure(closure))
        session.dag.addOperatorNode('merge', inputs, result)
        return result;
    }

    static <V> DataflowReadChannel<V> merge(final DataflowReadChannel source, final List<DataflowReadChannel<Object>> others, final Closure<V> closure) {
        final result = newChannelBy(source)
        final List<DataflowReadChannel<?>> inputs = new ArrayList<DataflowReadChannel<?>>(1 + others.size())
        inputs.add(source);
        inputs.addAll(others);
        newOperator(inputs, [result], new ChainWithClosure(closure));
        session.dag.addOperatorNode('merge', inputs, result)
        return result;
    }

}
