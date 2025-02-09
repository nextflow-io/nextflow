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

import static nextflow.extension.DataflowHelper.*
import static nextflow.splitter.SplitterFactory.*
import static nextflow.util.CheckHelper.*

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.NF
import nextflow.extension.op.ContextRunPerThread
import nextflow.extension.op.Op
import nextflow.script.ChannelOut
import nextflow.script.TokenBranchDef
import nextflow.script.TokenMultiMapDef
import nextflow.splitter.FastaSplitter
import nextflow.splitter.FastqSplitter
import nextflow.splitter.JsonSplitter
import nextflow.splitter.TextSplitter
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
/**
 * A set of operators inspired to RxJava extending the methods available on DataflowChannel
 * data structure
 *
 * See https://github.com/Netflix/RxJava/wiki/Observable
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class OperatorImpl {

    static final public OperatorImpl instance = new OperatorImpl()

    /**
     * Subscribe *onNext* event
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowReadChannel subscribe(final DataflowReadChannel source, final Closure closure) {
        new SubscribeOp()
            .withInput(source)
            .withOnNext(closure)
            .apply()
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
        new SubscribeOp()
            .withInput(source)
            .withEvents(events)
            .apply()
        return source
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them
     *
     * @param channel
     * @param closure
     * @return
     */
    DataflowWriteChannel map(final DataflowReadChannel source, final Closure closure) {
        return new MapOp()
            .withSource(source)
            .withMapper(closure)
            .apply()
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
        new FlatMapOp()
            .withSource(source)
            .withMapper(closure)
            .apply()
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

        ReduceOp .create()
            .withSource(source)
            .withAction(closure)
            .apply()
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

        ReduceOp .create()
            .withSource(source)
            .withSeed(seed)
            .withAction(closure)
            .apply()
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
        return new FilterOp()
            .withSource(source)
            .withCriteria(criteria)
            .apply()
    }

    DataflowWriteChannel filter(DataflowReadChannel source, final Closure<Boolean> closure) {
        return new FilterOp()
            .withSource(source)
            .withCriteria(closure)
            .apply()
    }

    DataflowWriteChannel until(DataflowReadChannel source, final Closure<Boolean> closure) {
        return new UntilOp(source,closure).apply()
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
        unique(source, Closure.IDENTITY)
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
    DataflowWriteChannel unique(final DataflowReadChannel source, final Closure comparator) {
        return new UniqueOp()
            .withSource(source)
            .withComparator(comparator)
            .apply()
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
        distinct(source, Closure.IDENTITY)
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    DataflowWriteChannel distinct(final DataflowReadChannel source, final Closure comparator) {
        new DistinctOp()
            .withSource(source)
            .withComparator(comparator)
            .apply()
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

        return first(source, { true })
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
    DataflowWriteChannel first( final DataflowReadChannel source, final Object criteria ) {
        new FirstOp()
            .withSource(source)
            .withCriteria(criteria)
            .apply()
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

        return new TakeOp(source,n).apply()
    }

    /**
     * The last operator creates a channel that only returns the last item emitted by the source channel
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the `last` item in the channel
     */
    DataflowWriteChannel last( final DataflowReadChannel source ) {
        new LastOp().withSource(source).apply()
    }

    DataflowWriteChannel collect(final DataflowReadChannel source, Closure action=null) {
        collect(source,Collections.emptyMap(),action)
    }

    DataflowWriteChannel collect(final DataflowReadChannel source, Map opts, Closure action=null) {
        return new CollectOp(source,action,opts).apply()
    }

    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    DataflowWriteChannel toList(final DataflowReadChannel source) {
        return new ToListOp(source).apply()
    }

    /**
     * Convert a {@code DataflowReadChannel} alias *channel* to a Java {@code List} sorting its content
     *
     * @param source The channel to be converted
     * @return A list holding all the items send over the channel
     */
    DataflowWriteChannel toSortedList(final DataflowReadChannel source, Closure closure = null) {
        return new ToListOp(source, closure ?: true).apply()
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param source
     * @param value
     * @return
     */
    DataflowWriteChannel count(final DataflowReadChannel source ) {
        return count0(source, null)
    }

    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param source
     * @param criteria
     * @return
     */
    DataflowWriteChannel count(final DataflowReadChannel source, final Object criteria ) {
        return count0(source, criteria)
    }

    private static DataflowVariable count0(DataflowReadChannel<?> source, Object criteria) {

        final target = new DataflowVariable()
        final discriminator = criteria != null ? new BooleanReturningMethodInvoker("isCase") : null

        final action = { current, item ->
            discriminator == null || discriminator.invoke(criteria, item) ? current+1 : current
        }

        ReduceOp .create()
            .withSource(source)
            .withTarget(target)
            .withSeed(0)
            .withAction(action)
            .apply()
        
        return target
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} returning the minimum value
     */
    DataflowWriteChannel min(final DataflowReadChannel source) {
        ReduceOp .create()
            .withSource(source)
            .withAction{ min, val -> val<min ? val : min }
            .apply()
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

        Closure action = null
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            action = (Closure){ min, item -> comparator.call(item) < comparator.call(min) ? item : min  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            action = (Closure){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        ReduceOp .create()
            .withSource(source)
            .withAction(action)
            .apply()
    }

    /**
     * The min operator waits until the source channel completes, and then emits the value that had the lowest value
     *
     * @param source The source channel
     * @param comparator The a {@code Comparator} object
     * @return A {@code DataflowVariable} returning the minimum value
     */
    DataflowWriteChannel min(final DataflowReadChannel source, Comparator comparator) {
        ReduceOp .create()
            .withSource(source)
            .withAction{ a, b -> comparator.compare(a,b)<0 ? a : b }
            .apply()
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    DataflowWriteChannel max(final DataflowReadChannel source) {
        ReduceOp .create()
            .withSource(source)
            .withAction { max, val -> val>max ? val : max }
            .apply()
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

        ReduceOp .create()
            .withSource(source)
            .withAction(action)
            .apply()
    }

    /**
     * The max operator waits until the source channel completes, and then emits the value that had the greatest value.
     *
     * @param source The source channel
     * @param comparator A {@code Comparator} object
     * @return A {@code DataflowVariable} emitting the maximum value
     */
    DataflowVariable max(final DataflowReadChannel source, Comparator comparator) {
        ReduceOp .create()
            .withSource(source)
            .withAction { a, b -> comparator.compare(a,b)>0 ? a : b }
            .apply()
    }

    /**
     * The sum operators crates a channel that emits the sum of all values emitted by the source channel to which is applied
     *
     * @param source  The source channel providing the values to sum
     * @param closure  A closure that given an entry returns the value to sum
     * @return A {@code DataflowVariable} emitting the final sum value
     */
    DataflowWriteChannel sum(final DataflowReadChannel source, Closure closure = null) {
        return MathOp.sum()
            .withSource(source)
            .withAction(closure)
            .apply()
    }


    DataflowWriteChannel mean(final DataflowReadChannel source, Closure closure = null) {
        return MathOp.mean()
            .withSource(source)
            .withAction(closure)
            .apply()
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
        new FlattenOp()
            .withSource(source)
            .apply()
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
        return new BufferOp(source)
                        .setParams(params)
                        .setCloseCriteria(closingCriteria)
                        .apply()
    }

    DataflowWriteChannel buffer( final DataflowReadChannel source, Object startingCriteria, Object closingCriteria ) {
        assert startingCriteria != null
        assert closingCriteria != null

        return new BufferOp(source)
                .setStartCriteria(startingCriteria)
                .setCloseCriteria(closingCriteria)
                .apply()
    }

    DataflowWriteChannel buffer( DataflowReadChannel source, Map<String,?> params ) {
        checkParams( 'buffer', params, 'size','skip','remainder' )
        return new BufferOp(source)
                        .setParams(params)
                        .apply()
    }

    DataflowWriteChannel collate( DataflowReadChannel source, int size, boolean keepRemainder = true ) {
        if( size <= 0 ) {
            throw new IllegalArgumentException("Illegal argument 'size' for operator 'collate' -- it must be greater than zero: $size")
        }

        new BufferOp(source)
            .setParams( size: size, remainder: keepRemainder )
            .apply()
    }

    DataflowWriteChannel collate( DataflowReadChannel source, int size, int step, boolean keepRemainder = true ) {
        new CollateOp()
            .withSource(source)
            .withSize(size)
            .withStep(step)
            .withRemainder(keepRemainder)
            .apply()
    }

    /**
     * Similar to https://github.com/Netflix/RxJava/wiki/Combining-Observables#merge
     *
     * @param source
     * @param others
     * @return
     */
    DataflowWriteChannel mix( DataflowReadChannel source, DataflowReadChannel[] others ) {
        if( others.size()==0 )
            throw new IllegalArgumentException("Operator 'mix' should have at least one right operand")

        return new MixOp(source,others).apply()
    }

    DataflowWriteChannel join( DataflowReadChannel left, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        // due to issue #582 the right channel cannot be provided in the join method signature
        // therefore the channel need to be added `'manually` to the inputs list
        // fixes #1346
        OpCall.current.get().inputs.add(right)
        return new JoinOp(left,right) .apply()
    }

    DataflowWriteChannel join( DataflowReadChannel left, Map opts, right ) {
        if( right==null ) throw new IllegalArgumentException("Operator `join` argument cannot be null")
        if( !(right instanceof DataflowReadChannel) ) throw new IllegalArgumentException("Invalid operator `join` argument [${right.getClass().getName()}] -- it must be a channel type")
        // due to issue #582 the right channel cannot be provided in the join method signature
        // therefore the channel need to be added `'manually` to the inputs list
        // fixes #1346
        OpCall.current.get().inputs.add(right)
        return new JoinOp(left,right,opts) .apply()
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
        return new CrossOp(source, other)
                    .setMapper(mapper)
                    .apply()
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
     * Implements a tap that create implicitly a new dataflow variable in the global script context.
     * For example:
     *
     * <pre>
     *     Channel.of(...)
     *            .tap { newChannelName }
     *            .map { ... }
     *  </pre>
     *
     * @param source The source dataflow variable
     * @param holder The closure defining the new variable name
     * @return The tap resulting dataflow channel
     */
    DataflowWriteChannel tap( final DataflowReadChannel source, final Closure holder ) {
        final tap = new TapOp(source, holder).apply()
        OpCall.current.get().outputs.addAll( tap.outputs )
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

        new SubscribeOp()
            .withInput(source)
            .withContext(new ContextRunPerThread())
            .withOnNext { DataflowProcessor dp, Object it -> Op.bind(dp,result,it); empty=false }
            .withOnComplete { DataflowProcessor dp ->
                    if(empty) {
                        final x = value instanceof Closure ? value.call() : value
                        Op.bind(dp,result,x)
                    }
                    if( !singleton )
                        result.bind(Channel.STOP)
                }
            .apply()

        return result
    }

    static private final Map PARAMS_VIEW = [newLine: Boolean]

    /**
     * Print out the channel content returning a new channel emitting the identical content as the original one
     *
     * @param source
     * @param closure
     * @return
     */
    DataflowWriteChannel view(final DataflowReadChannel source, Map opts, Closure closure = null) {
        assert source != null
        checkParams('view', opts, PARAMS_VIEW)
        final newLine = opts.newLine != false

        return new ViewOp()
            .withSource(source)
            .withNewLine(newLine)
            .withCode(closure)
            .apply()
    }

    DataflowWriteChannel view(final DataflowReadChannel source, Closure closure = null) {
        view(source, Collections.emptyMap(), closure)
    }

    // NO DAG
    DataflowWriteChannel merge(final DataflowReadChannel source, final DataflowReadChannel other, final Closure closure=null) {
        final result = CH.createBy(source)
        final inputs = [source, other]
        final action = closure ? new ChainWithClosure<>(closure) : new DefaultMergeClosure(inputs.size())
        final listener = stopErrorListener(source,result)
        new Op()
            .withInputs(inputs)
            .withOutput(result)
            .withListener(listener)
            .withCode(action)
            .apply()
        return result;
    }

    // NO DAG
    DataflowWriteChannel merge(final DataflowReadChannel source, final DataflowReadChannel... others) {
        new MergeOp(source,others as List).apply()
    }

    // NO DAG
    DataflowWriteChannel merge(final DataflowReadChannel source, final List<DataflowReadChannel> others, final Closure closure=null) {
        new MergeOp(source,others,closure).apply()
    }

    DataflowWriteChannel randomSample(DataflowReadChannel source, int n, Long seed = null) {
        if( source instanceof DataflowExpression )
            throw new IllegalArgumentException("Operator `randomSample` cannot be applied to a value channel")

        return new RandomSampleOp(source,n, seed).apply()
    }

    DataflowWriteChannel toInteger(final DataflowReadChannel source) {
        return new MapOp()
            .withSource(source)
            .withMapper { it -> it as Integer }
            .apply()
    }

    DataflowWriteChannel toLong(final DataflowReadChannel source) {
        return new MapOp()
            .withSource(source)
            .withMapper { it -> it as Long }
            .apply()
    }

    DataflowWriteChannel toFloat(final DataflowReadChannel source) {
        return new MapOp()
            .withSource(source)
            .withMapper { it -> it as Float }
            .apply()
    }

    DataflowWriteChannel toDouble(final DataflowReadChannel source) {
        return new MapOp()
            .withSource(source)
            .withMapper { it -> it as Double }
            .apply()
    }

    DataflowWriteChannel transpose( final DataflowReadChannel source, final Map params=null ) {
        return new TransposeOp(source,params).apply()
    }

    DataflowWriteChannel splitText(DataflowReadChannel source, Map opts=null) {
        return new SplitOp( source, 'splitText', opts ).apply()
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
        return new SplitOp( source, 'splitCsv', opts ).apply()
    }

    DataflowWriteChannel splitFasta(DataflowReadChannel source, Map opts=null) {
        return new SplitOp( source, 'splitFasta', opts ).apply()
    }

    DataflowWriteChannel splitFastq(DataflowReadChannel source, Map opts=null) {
        return new SplitOp( source, 'splitFastq', opts ).apply()
    }

    DataflowWriteChannel splitJson(DataflowReadChannel source, Map opts=null) {
        return new SplitOp( source, 'splitJson', opts ).apply()
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

    DataflowWriteChannel countJson(DataflowReadChannel source, Map opts=null) {
        final splitter = new JsonSplitter()
        final result = countOverChannel( source, splitter, opts )
        return result
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
            throw new IllegalArgumentException("Missing name with which to set the channel variable")

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

}
