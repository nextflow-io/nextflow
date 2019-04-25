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

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.splitter.FastaSplitter
import nextflow.splitter.FastqSplitter
import nextflow.splitter.TextSplitter
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
import static DataflowHelper.chainImpl
import static DataflowHelper.newOperator
import static DataflowHelper.reduceImpl
import static DataflowHelper.subscribeImpl
import static nextflow.extension.DataflowHelper.createOpParams
import static nextflow.extension.DataflowHelper.stopErrorListener
import static nextflow.splitter.SplitterFactory.countOverChannel
import static nextflow.util.CheckHelper.checkParams
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

    /**
     * Subscribe *onNext* event
     *
     * @param source
     * @param closure
     * @return
     */
    static final <V> DataflowReadChannel<V> subscribe(final DataflowReadChannel<V> source, final Closure<V> closure) {
        subscribeImpl( source, [onNext: closure] )
        NodeMarker.addOperatorNode('subscribe', source, null)
        return source
    }

    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    static final <V> DataflowReadChannel<V> subscribe(final DataflowReadChannel<V> source, final Map<String,Closure> events ) {
        subscribeImpl(source, events)
        NodeMarker.addOperatorNode('subscribe', source, null)
        return source
    }

    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Closure<V> closure) {

        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>(closure))
        NodeMarker.addOperatorNode('chain', source, target)

        return target;
    }


    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Map<String, Object> params, final Closure<V> closure) {

        final DataflowReadChannel<V> target = newChannelBy(source)
        chainImpl(source, target, params, closure)
        NodeMarker.addOperatorNode('chain', source, target)

        return target;
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them
     *
     * @param channel
     * @param closure
     * @return
     */
    static final <V> DataflowReadChannel<V> map(final DataflowReadChannel<?> source, final Closure<V> closure) {
        assert source != null
        assert closure

        def target = new MapOp(source, closure).apply()

        NodeMarker.addOperatorNode('map', source, target)
        return target;
    }


    /**
     * Transform the items emitted by a channel by applying a function to each of them and then flattens the results of that function.
     *
     * @param source The source channel
     * @param closure The closure mapping the values emitted by the source channel
     * @return The channel emitting the mapped values
     */
    static final <V> DataflowReadChannel<V> flatMap(final DataflowReadChannel<?> source, final Closure<V> closure=null) {
        assert source != null

        final target = new DataflowQueue()

        def listener = stopErrorListener(source,target)

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

        NodeMarker.addOperatorNode('flatMap', source, target)
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
    static final <V> DataflowReadChannel<V> reduce(final DataflowReadChannel<?> source, final Closure<V> closure) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException('Operator `reduce` cannot be applied to a value channel')

        final target = new DataflowVariable()
        reduceImpl( source, target, null, closure )
        NodeMarker.addOperatorNode('reduce', source, target)
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
    static final <V> DataflowReadChannel<V> reduce(final DataflowReadChannel<?> source, V seed, final Closure<V> closure) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException('Operator `reduce` cannot be applied to a value channel')

        final target = new DataflowVariable()
        reduceImpl( source, target, seed, closure )
        NodeMarker.addOperatorNode('reduce', source, target)
        return target
    }

    static final DataflowReadChannel collectFile( final DataflowReadChannel source, final Closure closure = null ) {
        final result = new CollectFileOp(source, null, closure).apply()
        NodeMarker.addOperatorNode('collectFile', source, result)
        return result
    }

    static final DataflowReadChannel collectFile( final DataflowReadChannel source, Map params, final Closure closure = null ) {
        def result = new CollectFileOp(source, params, closure).apply()
        NodeMarker.addOperatorNode('collectFile', source, result)
        return result
    }

    static final DataflowReadChannel groupTuple( final DataflowReadChannel source, final Map params ) {
        def result = new GroupTupleOp(params, source).apply()
        NodeMarker.addOperatorNode('groupTuple', source, result)
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
    static final <V> DataflowReadChannel<V> filter(final DataflowReadChannel<V> source, final Object criteria) {
        def discriminator = new BooleanReturningMethodInvoker("isCase");

        def target = newChannelBy(source)
        if( source instanceof DataflowExpression ) {
            source.whenBound {
                def result = it instanceof ControlMessage ? false : discriminator.invoke(criteria, (Object)it)
                target.bind( result ? it : Channel.STOP )
            }
        }
        else {
            newOperator(source, target, {
                def result = discriminator.invoke(criteria, (Object)it)
                if( result ) target.bind(it)
            })
        }

        NodeMarker.addOperatorNode('filter', source, target)
        return target
    }

    static <T> DataflowReadChannel<T> filter(DataflowReadChannel<T> source, final Closure<Boolean> closure) {
        def target = newChannelBy(source)
        if( source instanceof DataflowExpression ) {
            source.whenBound {
                def result = it instanceof ControlMessage ? false : DefaultTypeTransformation.castToBoolean(closure.call(it))
                target.bind( result ? it : Channel.STOP )
            }
        }
        else {
            newOperator(source, target, {
                def result = DefaultTypeTransformation.castToBoolean(closure.call(it))
                if( result ) target.bind(it)
            })
        }

        NodeMarker.addOperatorNode('filter', source, target)
        return target
    }

    static <T> DataflowReadChannel<T> until(DataflowReadChannel<T> source, final Closure<Boolean> closure) {
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
        NodeMarker.addOperatorNode('until', source, target)
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
    static final <V> DataflowReadChannel<V> unique(final DataflowReadChannel<V> source) {
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
    static final <V> DataflowReadChannel<V> unique(final DataflowReadChannel<V> source, Closure comparator ) {

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

        NodeMarker.addOperatorNode('unique', source, target)
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
    static final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> source ) {
        distinct(source) {it}
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    static final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> source, Closure<?> comparator ) {

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

        NodeMarker.addOperatorNode('distinct', source, target)
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
    static final <V> DataflowReadChannel<V> first( DataflowReadChannel<V> source ) {
        if( source instanceof DataflowExpression ) {
            def msg = "The operator `first` is useless when applied to a value channel which returns a single value by definition"
            def name = session?.binding?.getVariableName(source)
            if( name )
                msg += " -- check channel `$name`"
            log.warn msg
        }

        def target = new DataflowVariable<V>()
        source.whenBound { target.bind(it) }
        NodeMarker.addOperatorNode('first', source, target)
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
    static final <V> DataflowReadChannel<V> first( final DataflowReadChannel<V> source, Object criteria ) {

        def target = new DataflowVariable()
        def discriminator = new BooleanReturningMethodInvoker("isCase");

        if( source instanceof DataflowExpression ) {
            source.whenBound {
                def result = it instanceof ControlMessage ? false : discriminator.invoke(criteria, it)
                target.bind( result ? it : Channel.STOP )
            }
        }
        else {
            newOperator([source],[]) {
                if( discriminator.invoke(criteria, it) ) {
                    target.bind(it)
                    ((DataflowProcessor) getDelegate()).terminate()
                }
            }
        }

        NodeMarker.addOperatorNode('first', source, target)
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
    static final <V> DataflowReadChannel<V> take( final DataflowReadChannel<V> source, int n ) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException("Operator `take` cannot be applied to a value channel")

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

        NodeMarker.addOperatorNode('take', source, target)
        return target
    }

    /**
     * The last operator creates a channel that only returns the last item emitted by the source channel
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the `last` item in the channel
     */
    static final <V> DataflowReadChannel<V> last( final DataflowReadChannel<V> source ) {

        def target = new DataflowVariable()
        def V last = null
        subscribeImpl( source, [onNext: { last = it }, onComplete: {  target.bind(last) }] )
        NodeMarker.addOperatorNode('last', source, target)
        return target
    }

    static final <V> DataflowReadChannel<V> collect(final DataflowReadChannel<V> source, Closure action=null) {
        collect(source,Collections.emptyMap(),action)
    }

    static final <V> DataflowReadChannel<V> collect(final DataflowReadChannel<V> source, Map opts, Closure action=null) {
        final target = new CollectOp(source,action,opts).apply()
        NodeMarker.addOperatorNode('collect', source, target)
        return target
    }


    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static final <V> DataflowReadChannel<V> toList(final DataflowReadChannel<V> source) {
        final target = ToListOp.apply(source)
        NodeMarker.addOperatorNode('toList', source, target)
        return target
    }

    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List} sorting its content
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static final <V> DataflowReadChannel<V> toSortedList(final DataflowReadChannel<V> source, Closure closure = null) {
        final target = new ToListOp(source, closure ?: true).apply()
        NodeMarker.addOperatorNode('toSortedList', source, target)
        return target as DataflowVariable
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param source
     * @param value
     * @return
     */
    static final DataflowReadChannel<Number> count(final DataflowReadChannel<?> source ) {
        final target = count0(source, null)
        NodeMarker.addOperatorNode('count', source, target)
        return target
    }

    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param source
     * @param criteria
     * @return
     */
    static final DataflowReadChannel<Number> count(final DataflowReadChannel<?> source, final Object criteria ) {
        final target = count0(source, criteria)
        NodeMarker.addOperatorNode('count', source, target)
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
    static final DataflowReadChannel<Map> countBy(final DataflowReadChannel<?> source ) {
        countBy(source, { it })
    }

    /**
     * Sorts all collection members into groups determined by the supplied mapping closure and counts the group size
     *
     * @param source
     * @param criteria
     * @return
     */
    static final DataflowReadChannel<Map> countBy(final DataflowReadChannel<?> source, final Closure criteria ) {

        final target = new DataflowVariable()

        reduceImpl(source, target, [:]) { Map map, item ->
                def key = criteria.call(item)
                def value = map.containsKey(key) ? map.get(key)+1 : 1
                map.put(key, value)
                return map
        }

        NodeMarker.addOperatorNode('countBy', source, target)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} returning the minimum value
     */
    static final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> source) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { min, val -> val<min ? val : min }
        NodeMarker.addOperatorNode('min', source, target)
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
    static final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> source, Closure<V> comparator) {

        def action
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure<V>){ min, item -> comparator.call(item) < comparator.call(min) ? item : min  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure<V>){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        final target = new DataflowVariable()
        reduceImpl(source, target, null, action)
        NodeMarker.addOperatorNode('min', source, target)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @param comparator The a {@code Comparator} object
     * @return A {@code DataflowVariable} returning the minimum value
     */
    static final <V> DataflowReadChannel<V> min(final DataflowQueue<V> source, Comparator comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)<0 ? a : b }
        NodeMarker.addOperatorNode('min', source, target)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    static final <V> DataflowReadChannel<V> max(final DataflowQueue source) {
        final target = new DataflowVariable()
        reduceImpl(source,target, null) { max, val -> val>max ? val : max }
        NodeMarker.addOperatorNode('max', source, target)
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
    static final <V> DataflowReadChannel<V> max(final DataflowQueue<V> source, Closure comparator) {

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
        NodeMarker.addOperatorNode('max', source, target)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @param comparator A {@code Comparator} object
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    static final <V> DataflowVariable<V> max(final DataflowQueue<V> source, Comparator<V> comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)>0 ? a : b }
        NodeMarker.addOperatorNode('max', source, target)
        return target
    }

    /**
     * The sum operators crates a channel that emits the sum of all values emitted by the source channel to which is applied
     *
     * @param source  The source channel providing the values to sum
     * @param closure  A closure that given an entry returns the value to sum
     * @return A {@code DataflowVariable} emitting the final sum value
     */
    static final DataflowReadChannel sum(final DataflowQueue source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'sum', action: closure)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
        NodeMarker.addOperatorNode('sum', source, target)
        return target
    }


    static final DataflowReadChannel mean(final DataflowQueue source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'mean', action: closure, mean: true)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
        NodeMarker.addOperatorNode('mean', source, target)
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
    static final DataflowReadChannel<Map> groupBy(final DataflowReadChannel source, final params = null ) {
        log.warn "Operator `groupBy` is deprecated and it will be removed in a future release"
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

        NodeMarker.addOperatorNode('groupBy', source, target)
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
    static final void route( final DataflowReadChannel source, Map<?,DataflowWriteChannel> targets, Closure mapper = DEFAULT_MAPPING_CLOSURE ) {
        log.warn "Operator `route` is deprecated -- It will be removed in a future release"

        DataflowHelper.subscribeImpl(source,
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

        NodeMarker.addOperatorNode('route', source, targets.values())
    }

    static final DataflowReadChannel route( final DataflowReadChannel source, final Closure mapper = DEFAULT_MAPPING_CLOSURE ) {
        assert !(source instanceof DataflowExpression)
        log.warn "Operator `route` is deprecated and it will be removed in a future release"

        def allChannels = new ConcurrentHashMap()
        DataflowQueue target = new DataflowQueue()

        DataflowHelper.subscribeImpl(source,
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

        NodeMarker.addOperatorNode('route', source, target)
        return target
    }


    static final DataflowReadChannel spread( final DataflowReadChannel source, Object other ) {

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
            def left = [a]
            def right = (b instanceof List ? b : [b])
            [left, right]
                    .combinations()
                    .each{ Collection it -> proc.bindOutput(it.flatten())  }
        }

        def sources = new ArrayList(2)
        sources.add ( source )
        sources.add ( other instanceof DataflowChannel ? other : inputs  )
        NodeMarker.addOperatorNode('spread', sources, target)
        return target
    }

    static final DataflowReadChannel combine( DataflowReadChannel left, Object right ) {
        combine(left, null, right)
    }

    static final DataflowReadChannel combine( DataflowReadChannel left, Map params, Object right ) {
        checkParams('combine', params, [flat:Boolean, by: [List,Integer]])

        final op = new CombineOp(left,right)
        final sources = op.inputs
        if( params?.by != null ) op.pivot = params.by
        final target = op.apply()
        NodeMarker.addOperatorNode('combine', sources, target)
        return target
    }

    static final DataflowReadChannel flatten( final DataflowReadChannel source )  {

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

        NodeMarker.addOperatorNode('flatten', source, target)
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
    static final <V> DataflowReadChannel<V> buffer( final DataflowReadChannel<V> source, Map params=null, Object closingCriteria ) {

        def target = new BufferOp(source)
                        .setParams(params)
                        .setCloseCriteria(closingCriteria)
                        .apply()
        NodeMarker.addOperatorNode('buffer', source, target)
        return target
    }

    static final <V> DataflowReadChannel<V> buffer( final DataflowReadChannel<V> source, Object startingCriteria, Object closingCriteria ) {
        assert startingCriteria != null
        assert closingCriteria != null

        def target = new BufferOp(source)
                .setStartCriteria(startingCriteria)
                .setCloseCriteria(closingCriteria)
                .apply()

        NodeMarker.addOperatorNode('buffer', source, target)
        return target
    }

    static final <V> DataflowReadChannel<V> buffer( DataflowReadChannel<V> source, Map<String,?> params ) {
        checkParams( 'buffer', params, 'size','skip','remainder' )

        def target = new BufferOp(source)
                        .setParams(params)
                        .apply()

        NodeMarker.addOperatorNode('buffer', source, target)
        return target
    }


    static final <V> DataflowReadChannel<V> collate( DataflowReadChannel<V> source, int size, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        def target = new BufferOp(source)
                        .setParams( size: size, remainder: keepRemainder )
                        .apply()

        NodeMarker.addOperatorNode('collate', source, target)
        return target
    }

    static final <V> DataflowReadChannel<V> collate( DataflowReadChannel<V> source, int size, int step, boolean keepRemainder = true ) {
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

            Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
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

        NodeMarker.addOperatorNode('collate', source, target)
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

        DataflowHelper.subscribeImpl(source, handlers)
        others.each{ DataflowHelper.subscribeImpl(it, handlers) }

        def allSources = [source]
        allSources.addAll(others)

        NodeMarker.addOperatorNode('mix', allSources, target)
        return target
    }

    static final DataflowReadChannel join( DataflowReadChannel left, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        def target = new JoinOp(left,right) .apply()
        NodeMarker.addOperatorNode('join', [left, right], target)
        return target
    }

    static final DataflowReadChannel join( DataflowReadChannel left, Map opts, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        def target = new JoinOp(left,right,opts) .apply()
        NodeMarker.addOperatorNode('join', [left, right], target)
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

        NodeMarker.addOperatorNode('phase', [source, other], target)
        return target
    }

    static final DataflowReadChannel phase( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new PhaseOp(source,other)
                        .setMapper(mapper)
                        .apply()

        NodeMarker.addOperatorNode('phase', [source, other], target)
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

    static Closure DEFAULT_MAPPING_CLOSURE = { obj, int index=0 ->

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


    static <T> DataflowReadChannel cross( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new CrossOp(source, other)
                    .setMapper(mapper)
                    .apply()

        NodeMarker.addOperatorNode('cross', [source,other], target)
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

        NodeMarker.addOperatorNode('concat', allSources, target)
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
        NodeMarker.addOperatorNode('separate', source, outputs)
    }

    static void separate(final DataflowReadChannel source, final List<DataflowWriteChannel<?>> outputs) {
        new SeparateOp(source, outputs).apply()
        NodeMarker.addOperatorNode('separate', source, outputs)
    }

    static void separate(final DataflowReadChannel source, final List<DataflowWriteChannel<?>> outputs, final Closure<List<Object>> code) {
        new SeparateOp(source, outputs, code).apply()
        NodeMarker.addOperatorNode('separate', source, outputs)
    }

    static final List<DataflowReadChannel> separate( final DataflowReadChannel source, int n ) {
        def outputs = new SeparateOp(source, n).apply()
        NodeMarker.addOperatorNode('separate', source, outputs)
        return outputs
    }

    static final List<DataflowReadChannel> separate( final DataflowReadChannel source, int n, Closure mapper  ) {
        def outputs = new SeparateOp(source, n, mapper).apply()
        NodeMarker.addOperatorNode('separate', source, outputs)
        return outputs
    }


    static final void into( DataflowReadChannel source, final DataflowWriteChannel... targets ) {
        new IntoOp(source, targets as List<DataflowWriteChannel>).apply()
        NodeMarker.addOperatorNode('into', source, targets)
    }

    static final List<DataflowReadChannel> into( final DataflowReadChannel source, int n ) {
        def outputs = new IntoOp(source,n).apply().getOutputs()
        NodeMarker.addOperatorNode('into', source, outputs)
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
        NodeMarker.addOperatorNode('into', source, outputs)
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
        NodeMarker.addOperatorNode('tap', source, tap.outputs)
        return (DataflowReadChannel)tap.result
    }

    static DataflowReadChannel tap( final DataflowReadChannel source, final DataflowWriteChannel target ) {
        def tap = new TapOp(source, target).apply()
        NodeMarker.addOperatorNode('tap', source, tap.outputs)
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
    static DataflowReadChannel ifEmpty( DataflowReadChannel source, value ) {

        boolean empty = true
        final result = newChannelBy(source)
        final singleton = result instanceof DataflowExpression
        final next = { result.bind(it); empty=false }
        final complete = {
            if(empty)
                result.bind( value instanceof Closure ? value() : value )
            if( !singleton )
                result.bind(Channel.STOP)
        }

        subscribeImpl(source, [onNext: next, onComplete: complete])

        NodeMarker.addOperatorNode('ifEmpty', source, result)
        return result
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    static void print(final DataflowReadChannel<?> source, Closure closure = null) {
        final print0 = { def obj = closure ? closure.call(it) : it; session.printConsole(obj?.toString(),false) }
        subscribeImpl(source, [onNext: print0])
        NodeMarker.addOperatorNode('print', source, null)
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    static void println(final DataflowReadChannel<?> source, Closure closure = null) {
        final print0 = { def obj = closure ? closure.call(it) : it; session.printConsole(obj?.toString(),true) }
        subscribeImpl(source, [onNext: print0])
        NodeMarker.addOperatorNode('println', source, null)
    }


    static private final PARAMS_VIEW = [newLine: Boolean]

    /**
     * Print out the channel content retuning a new channel emitting the identical content as the original one
     *
     * @param source
     * @param closure
     * @return
     */
    static final <V> DataflowReadChannel<V> view(final DataflowReadChannel<?> source, Map opts, Closure closure = null) {
        assert source != null
        checkParams('view', opts, PARAMS_VIEW)
        final newLine = opts.newLine != false

        final target = newChannelBy(source);

        final apply = [

                onNext:
                        {
                            final obj = closure != null ? closure.call(it) : it
                            session.printConsole(obj?.toString(), newLine)
                            target.bind(it)
                        },

                onComplete: {
                    target.close()
                }
        ]

        subscribeImpl(source,apply)

        NodeMarker.addOperatorNode('view', source, target)
        return target;

    }

    static final <V> DataflowReadChannel<V> view(final DataflowReadChannel<?> source, Closure closure = null) {
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
        NodeMarker.addSourceNode('channel',target)
        return target
    }

    @Deprecated
    static DataflowBroadcast broadcast( DataflowReadChannel source ) {
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
        else if( source instanceof DataflowExpression ) {
            if( !source.isBound() )
                source.bind(Channel.STOP)
        }
        else {
            log.warn "Operator `close` cannot be applied to channels of type: ${source?.class?.simpleName}"
        }
        return source
    }

    static <T> void choice(final DataflowReadChannel source, final List<DataflowWriteChannel<T>> outputs, final Closure<Integer> code) {
        new ChoiceOp(source,outputs,code).apply()
        NodeMarker.addOperatorNode('choice', source, outputs)
    }

    static DataflowReadChannel merge(final DataflowReadChannel source, final DataflowReadChannel other, final Closure closure=null) {
        final result = newChannelBy(source)
        final inputs = [source, other]
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(inputs.size())
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, action)
        NodeMarker.addOperatorNode('merge', inputs, result)
        return result;
    }

    static DataflowReadChannel merge(final DataflowReadChannel source, final DataflowReadChannel... others) {
        final result = newChannelBy(source)
        final List<DataflowReadChannel> inputs = new ArrayList<DataflowReadChannel>(1 + others.size())
        inputs.add(source)
        inputs.addAll(others)
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, new DefaultMergeClosure(1 + others.size()))
        NodeMarker.addOperatorNode('merge', inputs, result)
        return result;
    }

    static DataflowReadChannel merge(final DataflowReadChannel source, final List<DataflowReadChannel> others, final Closure closure=null) {
        final result = newChannelBy(source)
        final List<DataflowReadChannel> inputs = new ArrayList<DataflowReadChannel>(1 + others.size())
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(1 + others.size())
        inputs.add(source)
        inputs.addAll(others)
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, action)
        NodeMarker.addOperatorNode('merge', inputs, result)
        return result;
    }

    static <V> DataflowReadChannel<V> randomSample(final DataflowReadChannel source, int n, Long seed = null) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException("Operator `randomSample` cannot be applied to a value channel")

        final result = new RandomSampleOp(source,n, seed).apply()
        NodeMarker.addOperatorNode('randomSample', source, result)
        return result;
    }

    static <V> DataflowReadChannel<V> toInteger(final DataflowReadChannel source) {
        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>({ it -> it as Integer }))
        NodeMarker.addOperatorNode('toInteger', source, target)
        return target;
    }

    static <V> DataflowReadChannel<V> toLong(final DataflowReadChannel source) {
        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>({ it -> it as Long }))
        NodeMarker.addOperatorNode('toLong', source, target)
        return target;
    }

    static <V> DataflowReadChannel<V> toFloat(final DataflowReadChannel source) {
        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>({ it -> it as Float }))
        NodeMarker.addOperatorNode('toFloat', source, target)
        return target;
    }

    static <V> DataflowReadChannel<V> toDouble(final DataflowReadChannel source) {
        final DataflowReadChannel<V> target = newChannelBy(source)
        newOperator(source, target, new ChainWithClosure<V>({ it -> it as Double }))
        NodeMarker.addOperatorNode('toDouble', source, target)
        return target;
    }

    static final DataflowReadChannel transpose( final DataflowReadChannel source, final Map params=null ) {
        def result = new TransposeOp(source,params).apply()
        NodeMarker.addOperatorNode('transpose', source, result)
        return result
    }

    static final DataflowReadChannel dump(final DataflowReadChannel source, Closure closure = null) {
        dump(source, Collections.emptyMap(), closure)
    }

    static final DataflowReadChannel dump(final DataflowReadChannel source, Map opts, Closure closure = null) {
        def op = new DumpOp(source, opts, closure)
        if( op.isEnabled() ) {
            def target = op.apply()
            NodeMarker.addOperatorNode('dump', source, target)
            return target;
        }
        else {
            return source
        }
    }


    static DataflowReadChannel splitText(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitText', opts ).apply()
        NodeMarker.addOperatorNode('splitText', source, result)
        return result
    }

    static DataflowReadChannel splitText(DataflowReadChannel source, Map opts=null, Closure action) {
        if( opts==null && action ) {
            opts = new HashMap<>(5)
        }
        opts.put('each', action)
        final result = new SplitOp( source, 'splitText', opts ).apply()
        NodeMarker.addOperatorNode('splitText', source, result)
        return result
    }

    static DataflowReadChannel splitCsv(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitCsv', opts ).apply()
        NodeMarker.addOperatorNode('splitCsv', source, result)
        return result
    }

    static DataflowReadChannel splitFasta(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitFasta', opts ).apply()
        NodeMarker.addOperatorNode('splitFasta', source, result)
        return result
    }

    static DataflowReadChannel splitFastq(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitFastq', opts ).apply()
        NodeMarker.addOperatorNode('splitFastq', source, result)
        return result
    }

    static DataflowReadChannel countLines(DataflowReadChannel source, Map opts=null) {
        final splitter = new TextSplitter()
        final result = countOverChannel( source, splitter, opts )
        NodeMarker.addOperatorNode('countLines', source, result)
        return result
    }

    static DataflowReadChannel countFasta(DataflowReadChannel source, Map opts=null) {
        final splitter = new FastaSplitter()
        final result = countOverChannel( source, splitter, opts )
        NodeMarker.addOperatorNode('countFasta', source, result)
        return result
    }

    static DataflowReadChannel countFastq(DataflowReadChannel source, Map opts=null) {
        final splitter = new FastqSplitter()
        final result = countOverChannel( source, splitter, opts )
        NodeMarker.addOperatorNode('countFastq', source, result)
        return result
    }


    static DataflowReadChannel countText(DataflowReadChannel source) {
        log.warn "Method `countText` has been deprecated -- Use `countLines` instead"
        countLines(source)
    }


}
