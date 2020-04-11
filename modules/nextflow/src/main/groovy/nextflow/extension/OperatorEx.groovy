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


import static nextflow.extension.DataflowHelper.*
import static nextflow.splitter.SplitterFactory.*
import static nextflow.util.CheckHelper.*

import java.lang.reflect.Modifier
import java.util.concurrent.atomic.AtomicInteger

import groovy.runtime.metaclass.DelegatingPlugin
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
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
import nextflow.NF
import nextflow.Session
import nextflow.script.ChannelOut
import nextflow.script.TokenBranchDef
import nextflow.script.TokenMultiMapDef
import nextflow.splitter.FastaSplitter
import nextflow.splitter.FastqSplitter
import nextflow.splitter.TextSplitter
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
class OperatorEx implements DelegatingPlugin {

    static final public OperatorEx instance = new OperatorEx()

    private static Session getSession() { Global.getSession() as Session }

    final public static Set<String> OPERATOR_NAMES

    static {
        OPERATOR_NAMES = getDeclaredExtensionMethods0()
        log.trace "Dataflow extension methods: ${OPERATOR_NAMES.sort().join(',')}"
    }

    @CompileStatic
    static private Set<String> getDeclaredExtensionMethods0() {
        def result = new HashSet<String>(30)
        def methods = OperatorEx.class.getDeclaredMethods()
        for( def handle : methods ) {
            if( !Modifier.isPublic(handle.getModifiers()) ) continue
            if( Modifier.isStatic(handle.getModifiers()) ) continue
            def params=handle.getParameterTypes()
            if( params.length>0 && isReadChannel(params[0]) )
                result.add(handle.name)
        }
        return result
    }

    @CompileStatic
    static boolean isReadChannel(Class clazz) {
        DataflowReadChannel.class.isAssignableFrom(clazz)
    }

    @CompileStatic
    boolean isExtensionMethod(Object obj, String name) {
        if( obj instanceof DataflowReadChannel || obj instanceof DataflowBroadcast || obj instanceof ChannelOut ) {
            return OPERATOR_NAMES.contains(name)
        }
        return false
    }


    @CompileStatic
    Object invokeExtensionMethod(Object channel, String method, Object[] args) {
        new OpCall(this,channel,method,args).call()
    }

    /**
     * Subscribe *onNext* event
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowReadChannel subscribe(final DataflowReadChannel source, final Closure closure) {
        subscribeImpl( source, [onNext: closure] )
        return source
    }

    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowReadChannel subscribe(final DataflowReadChannel source, final Map<String,Closure> events ) {
        subscribeImpl(source, events)
        return source
    }

    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowWriteChannel chain(final DataflowReadChannel<?> source, final Closure closure) {
        final target = CH.createBy(source)
        newOperator(source, target, new ChainWithClosure(closure))
        return target
    }

    /**
     * Chain operator, this is a synonym of {@code DataflowReadChannel.chainWith}
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowWriteChannel chain(final DataflowReadChannel<?> source, final Map<String, Object> params, final Closure closure) {
        final target = CH.createBy(source)
        chainImpl(source, target, params, closure)
        return target;
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them
     *
     * @param channel
     * @param closure
     * @return
     */
    DataflowWriteChannel map(final DataflowReadChannel source, final Closure closure) {
        assert source != null
        assert closure
        return new MapOp(source, closure).apply()
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them and then flattens the results of that function.
     *
     * @param source The source channel
     * @param closure The closure mapping the values emitted by the source channel
     * @return The channel emitting the mapped values
     */
    DataflowWriteChannel flatMap(final DataflowReadChannel<?> source, final Closure closure=null) {
        assert source != null

        final target = CH.create()
        final listener = stopErrorListener(source,target)

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
    DataflowWriteChannel reduce(final DataflowReadChannel source, final Closure closure) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException('Operator `reduce` cannot be applied to a value channel')

        final target = new DataflowVariable()
        reduceImpl( source, target, null, closure )
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
    DataflowWriteChannel reduce(final DataflowReadChannel<?> source, Object seed, final Closure closure) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException('Operator `reduce` cannot be applied to a value channel')

        final target = new DataflowVariable()
        reduceImpl( source, target, seed, closure )
        return target
    }

    DataflowWriteChannel collectFile( final DataflowReadChannel source, final Closure closure = null ) {
        final result = new CollectFileOp(source, null, closure).apply()
        return result
    }

    DataflowWriteChannel collectFile( final DataflowReadChannel source, Map params, final Closure closure = null ) {
        def result = new CollectFileOp(source, params, closure).apply()
        return result
    }

    DataflowWriteChannel groupTuple( final DataflowReadChannel source, final Map params=null ) {
        def result = new GroupTupleOp(params, source).apply()
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
    DataflowWriteChannel filter(final DataflowReadChannel source, final Object criteria) {
        def discriminator = new BooleanReturningMethodInvoker("isCase");

        def target = CH.createBy(source)
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

        return target
    }

    DataflowWriteChannel filter(DataflowReadChannel source, final Closure<Boolean> closure) {
        def target = CH.createBy(source)
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

        return target
    }

    DataflowWriteChannel until(DataflowReadChannel source, final Closure<Boolean> closure) {
        def target = CH.createBy(source)
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
    DataflowWriteChannel unique(final DataflowReadChannel source) {
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
    DataflowWriteChannel unique(final DataflowReadChannel source, Closure comparator ) {

        def history = [:]
        def target = CH.createBy(source)

        // when the operator stop clear the history map
        def events = new DataflowEventAdapter() {
            void afterStop(final DataflowProcessor processor) {
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
        }  as Closure

        // filter removing all duplicates
        chainImpl(source, target, [listeners: [events]], filter )

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
    DataflowWriteChannel distinct( final DataflowReadChannel source ) {
        distinct(source) {it}
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    DataflowWriteChannel distinct( final DataflowReadChannel source, Closure comparator ) {

        def previous = null
        final target = CH.createBy(source)
        Closure filter = { it ->

            def key = comparator.call(it)
            if( key == previous ) {
                return Channel.VOID
            }
            previous = key
            return it
        }

        chainImpl(source, target, [:], filter)

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
    DataflowWriteChannel first( DataflowReadChannel source ) {
        if( source instanceof DataflowExpression ) {
            def msg = "The operator `first` is useless when applied to a value channel which returns a single value by definition"
            def name = NF.lookupVariable(source)
            if( name )
                msg += " -- check channel `$name`"
            log.warn msg
        }

        def target = new DataflowVariable()
        source.whenBound { target.bind(it) }
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
    DataflowWriteChannel first( final DataflowReadChannel source, Object criteria ) {

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
    DataflowWriteChannel take( final DataflowReadChannel source, int n ) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException("Operator `take` cannot be applied to a value channel")

        def count = 0
        final target = CH.create()

        if( n==0 ) {
            target.bind(Channel.STOP)
            return target
        }

        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                if( ++count >= n ) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }
            }

            boolean onException(final DataflowProcessor processor, final Throwable e) {
                OperatorEx.log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        newOperator(
                inputs: [source],
                outputs: [target],
                listeners: (n > 0 ? [listener] : []),
                new ChainWithClosure(new CopyChannelsClosure()))

        return target
    }

    /**
     * The last operator creates a channel that only returns the last item emitted by the source channel
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the `last` item in the channel
     */
    DataflowWriteChannel last( final DataflowReadChannel source ) {

        def target = new DataflowVariable()
        def last = null
        subscribeImpl( source, [onNext: { last = it }, onComplete: {  target.bind(last) }] )
        return target
    }

    DataflowWriteChannel collect(final DataflowReadChannel source, Closure action=null) {
        collect(source,Collections.emptyMap(),action)
    }

    DataflowWriteChannel collect(final DataflowReadChannel source, Map opts, Closure action=null) {
        final target = new CollectOp(source,action,opts).apply()
        return target
    }


    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    DataflowWriteChannel toList(final DataflowReadChannel source) {
        final target = ToListOp.apply(source)
        return target
    }

    /**
     * Convert a {@code DataflowReadChannel} alias *channel* to a Java {@code List} sorting its content
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    DataflowWriteChannel toSortedList(final DataflowReadChannel source, Closure closure = null) {
        final target = new ToListOp(source, closure ?: true).apply()
        return target as DataflowVariable
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param source
     * @param value
     * @return
     */
    DataflowWriteChannel count(final DataflowReadChannel source ) {
        final target = count0(source, null)
        return target
    }

    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param source
     * @param criteria
     * @return
     */
    DataflowWriteChannel count(final DataflowReadChannel source, final Object criteria ) {
        final target = count0(source, criteria)
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
    @Deprecated
    DataflowWriteChannel<Map> countBy(final DataflowReadChannel source ) {
        countBy(source, { it })
    }

    /**
     * Sorts all collection members into groups determined by the supplied mapping closure and counts the group size
     *
     * @param source
     * @param criteria
     * @return
     */
    @Deprecated
    DataflowWriteChannel<Map> countBy(final DataflowReadChannel source, final Closure criteria ) {

        final target = new DataflowVariable()

        reduceImpl(source, target, [:]) { Map map, item ->
                def key = criteria.call(item)
                def value = map.containsKey(key) ? map.get(key)+1 : 1
                map.put(key, value)
                return map
        }

        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} returning the minimum value
     */
    DataflowWriteChannel min(final DataflowReadChannel source) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { min, val -> val<min ? val : min }
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
    DataflowWriteChannel min(final DataflowReadChannel source, Closure comparator) {

        def action
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure){ min, item -> comparator.call(item) < comparator.call(min) ? item : min  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        final target = new DataflowVariable()
        reduceImpl(source, target, null, action)
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @param comparator The a {@code Comparator} object
     * @return A {@code DataflowVariable} returning the minimum value
     */
    DataflowWriteChannel min(final DataflowReadChannel source, Comparator comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)<0 ? a : b }
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    DataflowWriteChannel max(final DataflowReadChannel source) {
        final target = new DataflowVariable()
        reduceImpl(source,target, null) { max, val -> val>max ? val : max }
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
    DataflowWriteChannel max(final DataflowReadChannel source, Closure comparator) {

        def action
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure){ max, item -> comparator.call(item) > comparator.call(max) ? item : max  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure){ a, b ->  comparator.call(a,b)>0 ? a : b  }
        }
        else {
            throw new IllegalArgumentException("Comparator closure can accept at most 2 arguments")
        }

        final target = new DataflowVariable()
        reduceImpl(source, target, null, action)
        return target
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @param comparator A {@code Comparator} object
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    DataflowVariable max(final DataflowReadChannel source, Comparator comparator) {
        final target = new DataflowVariable()
        reduceImpl(source, target, null) { a, b -> comparator.compare(a,b)>0 ? a : b }
        return target
    }

    /**
     * The sum operators crates a channel that emits the sum of all values emitted by the source channel to which is applied
     *
     * @param source  The source channel providing the values to sum
     * @param closure  A closure that given an entry returns the value to sum
     * @return A {@code DataflowVariable} emitting the final sum value
     */
    DataflowWriteChannel sum(final DataflowReadChannel source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'sum', action: closure)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
        return target
    }


    DataflowWriteChannel mean(final DataflowReadChannel source, Closure closure = null) {

        def target = new DataflowVariable()
        def aggregate = new Aggregate(name: 'mean', action: closure, mean: true)
        subscribeImpl(source, [onNext: aggregate.&process, onComplete: { target.bind( aggregate.result ) }])
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
    @Deprecated
    DataflowWriteChannel<Map> groupBy(final DataflowReadChannel source, final params = null ) {
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
        reduceImpl(source, target, [:]) { Map map, item ->
            def key = len == 2 ? mapper.call(item,index) : mapper.call(item)
            def list = map.get(key)
            list = list ? list << item : [item]
            map.put(key, list)
            return map
        }

        return target
    }

    @Deprecated
    DataflowWriteChannel spread( final DataflowReadChannel source, Object other ) {

        final target = CH.create()

        def inputs
        switch(other) {
            case DataflowExpression:
                inputs = other
                break
            case DataflowReadChannel:
                inputs = ToListOp.apply((DataflowReadChannel)other);
                break
            case Collection:
                inputs = Channel.value(other)
                OpCall.current.get().inputs.add(inputs)
                break
            case (Object[]):
                inputs = Channel.value(other as List)
                OpCall.current.get().inputs.add(inputs)
                break
            default:
                throw new IllegalArgumentException("Not a valid argument for 'spread' operator [${other?.class?.simpleName}]: ${other} -- Use a Collection or a channel instead. ")
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
            boolean onException(final DataflowProcessor processor, final Throwable e) {
                OperatorEx.log.error("@unknown", e)
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

        return target
    }

    DataflowWriteChannel combine( DataflowReadChannel left, Object right ) {
        combine(left, null, right)
    }

    DataflowWriteChannel combine( DataflowReadChannel left, Map params, Object right ) {
        checkParams('combine', params, [flat:Boolean, by: [List,Integer]])

        final op = new CombineOp(left,right)
        OpCall.current.get().inputs.addAll(op.inputs)
        if( params?.by != null ) op.pivot = params.by
        final target = op.apply()
        return target
    }

    DataflowWriteChannel flatten( final DataflowReadChannel source )  {

        final listeners = []
        final target = CH.create()

        if( source instanceof DataflowExpression ) {
            listeners << new DataflowEventAdapter() {
                @Override
                void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }

                boolean onException(final DataflowProcessor processor, final Throwable e) {
                    OperatorEx.log.error("@unknown", e)
                    session.abort(e)
                    return true;
                }
            }
        }


        newOperator(inputs: [source], outputs: [target], listeners: listeners) {  item ->

            def proc = ((DataflowProcessor) getDelegate())
            switch( item ) {
                case Collection:
                    ((Collection)item).flatten().each { value -> proc.bindOutput(value) }
                    break

                case (Object[]):
                    ((Collection)item).flatten().each { value -> proc.bindOutput(value) }
                    break

                case Channel.VOID:
                    break

                default:
                    proc.bindOutput(item)
            }
        }

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
    DataflowWriteChannel buffer( final DataflowReadChannel source, Map params=null, Object closingCriteria ) {

        def target = new BufferOp(source)
                        .setParams(params)
                        .setCloseCriteria(closingCriteria)
                        .apply()
        return target
    }

    DataflowWriteChannel buffer( final DataflowReadChannel source, Object startingCriteria, Object closingCriteria ) {
        assert startingCriteria != null
        assert closingCriteria != null

        def target = new BufferOp(source)
                .setStartCriteria(startingCriteria)
                .setCloseCriteria(closingCriteria)
                .apply()

        return target
    }

    DataflowWriteChannel buffer( DataflowReadChannel source, Map<String,?> params ) {
        checkParams( 'buffer', params, 'size','skip','remainder' )

        def target = new BufferOp(source)
                        .setParams(params)
                        .apply()

        return target
    }


    DataflowWriteChannel collate( DataflowReadChannel source, int size, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        def target = new BufferOp(source)
                        .setParams( size: size, remainder: keepRemainder )
                        .apply()

        return target
    }

    DataflowWriteChannel collate( DataflowReadChannel source, int size, int step, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        if( step <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'step' for operator 'collate' -- it must be greater than zero: $step")
        }

        // the result queue
        final target = CH.create()

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
                OperatorEx.log.error("@unknown", e)
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

        return target
    }


    /**
     * Similar to https://github.com/Netflix/RxJava/wiki/Combining-Observables#merge
     *
     * @param source
     * @param others
     * @return
     */
    DataflowWriteChannel mix( DataflowReadChannel source, DataflowReadChannel[] others ) {
        assert others.size()>0

        def target = CH.create()
        def count = new AtomicInteger( others.size()+1 )
        def handlers = [
                onNext: { target << it },
                onComplete: { if(count.decrementAndGet()==0) { target << Channel.STOP } }
        ]

        subscribeImpl(source, handlers)
        for( def it : others ) {
            subscribeImpl(it, handlers)
        }

        def allSources = [source]
        allSources.addAll(others)

        return target
    }

    DataflowWriteChannel join( DataflowReadChannel left, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        // due to issue #582 the right channel cannot be provided in the join method signature
        // therefore the channel need to be added `'manually` to the inputs list
        // fixes #1346
        OpCall.current.get().inputs.add(right)
        def target = new JoinOp(left,right) .apply()
        return target
    }

    DataflowWriteChannel join( DataflowReadChannel left, Map opts, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        // due to issue #582 the right channel cannot be provided in the join method signature
        // therefore the channel need to be added `'manually` to the inputs list
        // fixes #1346
        OpCall.current.get().inputs.add(right)
        def target = new JoinOp(left,right,opts) .apply()
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
    DataflowWriteChannel phase( DataflowReadChannel source, Map opts, DataflowReadChannel other, Closure mapper = null ) {

        def target = new PhaseOp(source,other)
                        .setMapper(mapper)
                        .setOpts(opts)
                        .apply()

        return target
    }

    DataflowWriteChannel phase( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new PhaseOp(source,other)
                        .setMapper(mapper)
                        .apply()

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

    @PackageScope
    static public Closure DEFAULT_MAPPING_CLOSURE = { obj, int index=0 ->

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


    DataflowWriteChannel cross( DataflowReadChannel source, DataflowReadChannel other, Closure mapper = null ) {

        def target = new CrossOp(source, other)
                    .setMapper(mapper)
                    .apply()

        return target
    }


    /**
     * Creates a channel that emits the items in same order as they are emitted by two or more channel
     *
     * @param source
     * @param others
     * @return
     */
    DataflowWriteChannel concat( DataflowReadChannel source, DataflowReadChannel... others ) {

        def target = new ConcatOp(source, others).apply()

        def allSources = [source]
        if(others) allSources.addAll(others)

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
    @DeprecatedDsl2
    void separate( DataflowReadChannel source, final DataflowWriteChannel... outputs ) {
        new SeparateOp(source, outputs as List<DataflowWriteChannel>).apply()
    }

    @DeprecatedDsl2
    void separate(final DataflowReadChannel source, final List<DataflowWriteChannel> outputs) {
        new SeparateOp(source, outputs).apply()
    }

    @DeprecatedDsl2
    void separate(final DataflowReadChannel source, final List<DataflowWriteChannel> outputs, final Closure<List> code) {
        new SeparateOp(source, outputs, code).apply()
    }

    @DeprecatedDsl2
    List<DataflowReadChannel> separate( final DataflowReadChannel source, int n ) {
        def outputs = new SeparateOp(source, n).apply()
        return outputs
    }

    @DeprecatedDsl2
    List<DataflowReadChannel> separate( final DataflowReadChannel source, int n, Closure mapper  ) {
        def outputs = new SeparateOp(source, n, mapper).apply()
        return outputs
    }

    @DeprecatedDsl2
    void into( DataflowReadChannel source, final DataflowWriteChannel... targets ) {
        new IntoOp(source, targets as List<DataflowWriteChannel>).apply()
    }

    @DeprecatedDsl2
    List<DataflowReadChannel> into( final DataflowReadChannel source, int n ) {
        def outputs = new IntoOp(source,n).apply().getOutputs()
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
    @DeprecatedDsl2
    void into( DataflowReadChannel source, Closure holder ) {
        def outputs = new IntoOp(source,holder).apply().getOutputs()
        OpCall.current.get().outputs.addAll(outputs)
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
    DataflowWriteChannel tap( final DataflowReadChannel source, final Closure holder ) {
        def tap = new TapOp(source, holder).apply()
        OpCall.current.get().outputs.addAll( tap.outputs )
        return tap.result
    }

    DataflowWriteChannel tap( final DataflowReadChannel source, final DataflowWriteChannel target ) {
        def tap = new TapOp(source, target).apply()
        return tap.result
    }


    /**
     * Empty the specified value only if the source channel to which is applied is empty i.e. do not emit
     * any value.
     *
     * @param source The channel to which the operator is applied
     * @param value The value to emit when the source channel is empty. If a closure is used the the value returned by its invocation is used.
     * @return The resulting channel emitting the source items or the default value when the channel is empty
     */
    DataflowWriteChannel ifEmpty( DataflowReadChannel source, value ) {

        boolean empty = true
        final result = CH.createBy(source)
        final singleton = result instanceof DataflowExpression
        final next = { result.bind(it); empty=false }
        final complete = {
            if(empty)
                result.bind( value instanceof Closure ? value() : value )
            if( !singleton )
                result.bind(Channel.STOP)
        }

        subscribeImpl(source, [onNext: next, onComplete: complete])

        return result
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    @DeprecatedDsl2(message='Operator `print` is deprecated -- Use `view` instead')
    void print(final DataflowReadChannel<?> source, Closure closure = null) {
        final print0 = { def obj = closure ? closure.call(it) : it; session.printConsole(obj?.toString(),false) }
        subscribeImpl(source, [onNext: print0])
    }

    /**
     * Print the channel content to the console standard output
     * @param source
     * @param closure
     */
    @DeprecatedDsl2(message='Operator `println` is deprecated -- Use `view` instead')
    void println(final DataflowReadChannel<?> source, Closure closure = null) {
        final print0 = { def obj = closure ? closure.call(it) : it; session.printConsole(obj?.toString(),true) }
        subscribeImpl(source, [onNext: print0])
    }


    static private final Map PARAMS_VIEW = [newLine: Boolean]

    /**
     * Print out the channel content retuning a new channel emitting the identical content as the original one
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowWriteChannel view(final DataflowReadChannel source, Map opts, Closure closure = null) {
        assert source != null
        checkParams('view', opts, PARAMS_VIEW)
        final newLine = opts.newLine != false

        final target = CH.createBy(source);
        final apply = new HashMap<String,Closure>(2)

        apply.onNext  = {
            final obj = closure != null ? closure.call(it) : it
            session.printConsole(obj?.toString(), newLine)
            target.bind(it)
        }

        apply. onComplete = { CH.close0(target) }

        subscribeImpl(source,apply)
        return target
    }

    DataflowWriteChannel view(final DataflowReadChannel source, Closure closure = null) {
        view(source, Collections.emptyMap(), closure)
    }

    void choice(final DataflowReadChannel source, final List<DataflowWriteChannel> outputs, final Closure<Integer> code) {
        new ChoiceOp(source,outputs,code).apply()
    }

    // NO DAG
    @DeprecatedDsl2
    DataflowWriteChannel merge(final DataflowReadChannel source, final DataflowReadChannel other, final Closure closure=null) {
        final result = CH.createBy(source)
        final inputs = [source, other]
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(inputs.size())
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, action)
        return result;
    }

    // NO DAG
    @DeprecatedDsl2
    DataflowWriteChannel merge(final DataflowReadChannel source, final DataflowReadChannel... others) {
        final result = CH.createBy(source)
        final List<DataflowReadChannel> inputs = new ArrayList<DataflowReadChannel>(1 + others.size())
        inputs.add(source)
        inputs.addAll(others)
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, new DefaultMergeClosure(1 + others.size()))
        return result;
    }

    // NO DAG
    @DeprecatedDsl2
    DataflowWriteChannel merge(final DataflowReadChannel source, final List<DataflowReadChannel> others, final Closure closure=null) {
        final result = CH.createBy(source)
        final List<DataflowReadChannel> inputs = new ArrayList<DataflowReadChannel>(1 + others.size())
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(1 + others.size())
        inputs.add(source)
        inputs.addAll(others)
        final listener = stopErrorListener(source,result)
        final params = createOpParams(inputs, result, listener)
        newOperator(params, action)
        return result;
    }

    DataflowWriteChannel randomSample(DataflowReadChannel source, int n, Long seed = null) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException("Operator `randomSample` cannot be applied to a value channel")

        final result = new RandomSampleOp(source,n, seed).apply()
        return result
    }

    DataflowWriteChannel toInteger(final DataflowReadChannel source) {
        final target = CH.createBy(source)
        newOperator(source, target, new ChainWithClosure({ it -> it as Integer }))
        return target;
    }

    DataflowWriteChannel toLong(final DataflowReadChannel source) {
        final target = CH.createBy(source)
        newOperator(source, target, new ChainWithClosure({ it -> it as Long }))
        return target;
    }

    DataflowWriteChannel toFloat(final DataflowReadChannel source) {
        final target = CH.createBy(source)
        newOperator(source, target, new ChainWithClosure({ it -> it as Float }))
        return target;
    }

    DataflowWriteChannel toDouble(final DataflowReadChannel source) {
        final target = CH.createBy(source)
        newOperator(source, target, new ChainWithClosure({ it -> it as Double }))
        return target;
    }

    DataflowWriteChannel transpose( final DataflowReadChannel source, final Map params=null ) {
        def result = new TransposeOp(source,params).apply()
        return result
    }

    DataflowWriteChannel splitText(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitText', opts ).apply()
        return result
    }

    DataflowWriteChannel splitText(DataflowReadChannel source, Map opts=null, Closure action) {
        if( opts==null && action ) {
            opts = new HashMap<>(5)
        }
        opts.put('each', action)
        final result = new SplitOp( source, 'splitText', opts ).apply()
        return result
    }

    DataflowWriteChannel splitCsv(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitCsv', opts ).apply()
        return result
    }

    DataflowWriteChannel splitFasta(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitFasta', opts ).apply()
        return result
    }

    DataflowWriteChannel splitFastq(DataflowReadChannel source, Map opts=null) {
        final result = new SplitOp( source, 'splitFastq', opts ).apply()
        return result
    }

    DataflowWriteChannel countLines(DataflowReadChannel source, Map opts=null) {
        final splitter = new TextSplitter()
        final result = countOverChannel( source, splitter, opts )
        return result
    }

    DataflowWriteChannel countFasta(DataflowReadChannel source, Map opts=null) {
        final splitter = new FastaSplitter()
        final result = countOverChannel( source, splitter, opts )
        return result
    }

    DataflowWriteChannel countFastq(DataflowReadChannel source, Map opts=null) {
        final splitter = new FastqSplitter()
        final result = countOverChannel( source, splitter, opts )
        return result
    }

    @Deprecated
    DataflowWriteChannel countText(DataflowReadChannel source) {
        countLines(source)
    }

    /**
     * Implement a `set` operator e.g.
     * <pre>
     *     SomeProcess.out | set { channelName }
     * </pre>
     *
     * @param source The channel instance to be bound in the context
     * @param holder A closure defining a variable identifier
     */

    void set(DataflowReadChannel source, Closure holder) {
        set0(source, holder)
    }

    void set(DataflowBroadcast source, Closure holder) {
        set0(source, holder)
    }

    void set(ChannelOut source, Closure holder) {
        set0(source, holder)
    }

    private void set0(source, Closure holder) {
        final name = CaptureProperties.capture(holder)
        if( !name )
            throw new IllegalArgumentException("Missing name to which set the channel variable")

        if( name.size()>1 )
            throw new IllegalArgumentException("Operation `set` does not allow more than one target name")

        NF.binding.setVariable(name[0], source)
        // do not add this nod in the DAG because it's not a real operator
        // since it's not transforming the channel
        OpCall.current.get().ignoreDagNode = true
    }

    /**
     * Implements `branch` operator
     *
     * @param source
     * @param action
     * @return
     */
    ChannelOut branch(DataflowReadChannel source, Closure<TokenBranchDef> action) {
        new BranchOp(source, action)
                .apply()
                .getOutput()
    }

    ChannelOut multiMap(DataflowReadChannel source, Closure<TokenMultiMapDef> action) {
        new MultiMapOp(source, action)
                .apply()
                .getOutput()
    }

    @Deprecated
    ChannelOut fork(DataflowReadChannel source, Closure<TokenMultiMapDef> action) {
        log.warn "Operator `fork` has been renamed to `multiMap`"
        multiMap(source, action)
    }
}
