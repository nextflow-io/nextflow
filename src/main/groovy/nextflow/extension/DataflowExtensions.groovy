package nextflow.extension
import static java.util.Arrays.asList

import java.util.concurrent.atomic.AtomicInteger

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
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.SeparationClosure
import nextflow.Channel
import org.codehaus.groovy.runtime.NullObject
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
class DataflowExtensions {


    /**
     * Create a dataflow object by the type of the specified source argument
     *
     * @param source
     * @return
     */
    static private final <V> DataflowChannel<V> newChannelBy(DataflowReadChannel<?> source) {

        switch( source ) {
            case DataflowExpression:
                return new DataflowVariable<V>()

            case DataflowQueue:
                return new DataflowQueue<V>()

            default:
                throw new IllegalArgumentException()
        }

    }

    static private VALID_HANDLERS = [ 'onNext', 'onComplete', 'onError' ]

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
        subscribe( source, [onNext: closure] )
    }


    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param source
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> subscribe(final DataflowReadChannel<V> source, final Map<String,Closure> events ) {
        checkSubscribeHandlers(events)

        def stopOnFirst = source instanceof DataflowExpression

        def listeners = []

        if( events.onNext ) {
            listeners << new DataflowEventAdapter() {
                @Override
                public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                    events.onNext.call(messages[0])
                    if( stopOnFirst ) {
                        processor.terminate()
                    }
                }
            }
        }

        if( events.onComplete ) {
            listeners << new DataflowEventAdapter() {
                @Override
                public void afterStop(final DataflowProcessor processor) {
                    events.onComplete.call(processor)
                }
            }
        }

        if( events.onException ) {
            listeners << new DataflowEventAdapter() {
                @Override
                public boolean onException(final DataflowProcessor processor, final Throwable e) {
                    events.onException.call(e)
                }
            }
        }


        final DataflowReadChannel<V> target = newChannelBy(source)
        final Map<String, Object> parameters = new HashMap<String, Object>();
        parameters.put("inputs", [source])
        parameters.put("outputs", [target])
        parameters.put('listeners', listeners)

        Dataflow.operator(parameters, {it});

        // forwards all
        return source

    }


    public static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Closure<V> closure) {
        final DataflowReadChannel<V> target = newChannelBy(source)
        Dataflow.operator(source, target, new ChainWithClosure<V>(closure))
        return target;
    }


    public static <V> DataflowReadChannel<V> chain(final DataflowReadChannel<?> source, final Map<String, Object> params, final Closure<V> closure) {

        final DataflowReadChannel<V> target = newChannelBy(source)
        final Map<String, Object> parameters = new HashMap<String, Object>(params)
        parameters.put("inputs", asList(source))
        parameters.put("outputs", asList(target))

        Dataflow.operator(parameters, new ChainWithClosure<V>(closure))

        return target;
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

        DataflowReadChannel<V> target = newChannelBy(source);
        Dataflow.operator(source, target, new ChainWithClosure<V>(closure));
        return target;

    }

    static public final <V> DataflowReadChannel<V> mapWithIndex( DataflowReadChannel<?> source, final Closure<V> closure ) {
        assert source != null
        assert closure

        int index = 0
        DataflowReadChannel<V> target = newChannelBy(source)

        Dataflow.operator(source, target) {

            def value = closure.call(it, index++)
            if (value != NullObject.getNullObject()) {
                ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(value);
            }
        }

        return target
    }

    /**
     * Transform the items emitted by a channel by applying a function to each of them and then flattens the results of that function.
     *
     * @param source The source channel
     * @param closure The closure mapping the values emitted by the source channel
     * @return The channel emitting the mapped values
     */
    static public final <V> DataflowReadChannel<V> mapMany(final DataflowReadChannel<?> source, final Closure<V> closure) {
        assert source != null
        assert closure

        def target = newChannelBy(source)

        Dataflow.operator(source, target) {

            def value = closure.call(it);
            switch( value ) {
                case Iterable:
                    ((Iterable)value) .each { entry -> bindOutput(entry) }
                    break

                case Map:
                    ((Map)value) .each { entry -> bindOutput(entry) }
                    break

                default:
                    bindOutput(value)
            }
        }

        return target
    }

    static private UNDEF = new Object()

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
        reduceImpl( source, null, closure )
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
        reduceImpl( source, seed, closure )
    }

    /**
     *
     *
     * @param channel
     * @param seed
     * @param closure
     * @return
     */
    static private <V> DataflowReadChannel<V> reduceImpl(final DataflowReadChannel<?> channel, def seed, final Closure<V> closure) {

        // the dataflow variable to return the final aggregation value
        def result = new DataflowVariable()

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
                DataflowExtensions.log.error('Dataflow reduce() exception', e)
                return true;
            }
        }


        channel.chain(listeners: [listener], {true})
        return result
    }

    /**
     * Iterates over the collection of items and returns each item that matches the given filter
     * by calling the {@code Object#isCase}method used by switch statements.
     *
     * This method can be used with different kinds of filters like regular expressions, classes, ranges etc. Example:
     *
     * def list = ['a', 'b', 'aa', 'bc', 3, 4.5]
     * assert list.grep( ~/a+/ )  == ['a', 'aa']
     * assert list.grep( ~/../ )  == ['aa', 'bc']
     * assert list.grep( Number ) == [ 3, 4.5 ]
     * assert list.grep{ it.toString().size() == 1 } == [ 'a', 'b', 3 ]
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final <V> DataflowReadChannel<V> grep(final DataflowReadChannel<V> channel, final Object criteria) {
        def discriminator = new BooleanReturningMethodInvoker("isCase");
        return channel.filter { Object it -> discriminator.invoke(criteria, it) }
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
        return source.chain(listeners: [events], filter)

    }

    /**
     * Suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     *
     * @return
     */
    static public final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> channel ) {
        distinct(channel) {it}
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    static public final <V> DataflowReadChannel<V> distinct( final DataflowReadChannel<V> channel, Closure<?> comparator ) {

        def previous = null

        return channel.chain { it ->

            def key = comparator.call(it)
            if( key == previous ) {
                return NullObject.getNullObject()
            }
            else {
                previous = key
                return it
            }
        }

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

        def discriminator = new BooleanReturningMethodInvoker("isCase");

        def result = new DataflowVariable()
        Dataflow.operator([source],[]) {

            if( discriminator.invoke(criteria, it) ) {
                result.bind(it)
                terminate()
            }
        }

        return result
    }

    /**
     *
     * emit only the first n items emitted by an Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#take
     *
     * @param source
     * @return
     */
    static public final <V> DataflowReadChannel<V> take( final DataflowReadChannel<V> source, int n ) {
        assert !(source instanceof DataflowExpression)

        def count = 0
        def target = new DataflowQueue<V>()
        Dataflow.operator([source],[]) {

            if( count++ < n ) {
                target << it
                return
            }
            target << Channel.STOP
            terminate()
        }

        return target
    }


    static public final <V> DataflowReadChannel<V> last( final DataflowReadChannel<V> source  ) {

        def target = new DataflowVariable()
        def V last = null
        source.subscribe( onNext: { last = it }, onComplete: {  target.bind(last) } )
        return target

    }

//
//
//    /**
//     * Emit items from a source Observable, but issue an exception if no item is emitted in a specified timespan
//     *
//     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#timeout
//     *
//     * @param channel
//     */
//    static public final <V> DataflowReadChannel<V> timeout( final DataflowQueue channel, String duration ) {
//        def millis = duration.isLong() ? duration.toLong() : Duration.create(duration).toMillis()
//        timeout( channel, millis )
//    }
//
//    static public final <V> DataflowReadChannel<V> timeout( final DataflowQueue channel, long duration, TimeUnit unit = TimeUnit.MILLISECONDS ) {
//
//    }

    /**
     * See https://github.com/Netflix/RxJava/wiki/Observable-Utility-Operators#parallel
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> parallel( final DataflowQueue<V> channel, Closure closure ) {

        def nProcs = Math.max( 2, Runtime.getRuntime().availableProcessors()-1)

        channel.filter( maxForks: nProcs ) { true }

    }

    /**
     * See  https://github.com/Netflix/RxJava/wiki/Observable-Utility-Operators#finallydo
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> doFinally( final DataflowReadChannel<V> channel, Closure<V> closure ) {
        channel.subscribe( onComplete: closure )
    }

    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param channel The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static public final <V> DataflowReadChannel<V> toList(final DataflowReadChannel<V> channel) {
        return reduce(channel, []) { list, item -> list << item }
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param channel
     * @param value
     * @return
     */
    static public final DataflowReadChannel<Number> count(final DataflowReadChannel<?> channel ) {
        reduce(channel, 0) { current, item -> current+1 }
    }


    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param source
     * @param criteria
     * @return
     */
    static public final DataflowReadChannel<Number> count(final DataflowReadChannel<?> source, final Object criteria ) {

        def discriminator = new BooleanReturningMethodInvoker("isCase");

        reduce(source, 0) { current, item ->
            discriminator.invoke(criteria, item) ? current+1 : current
        }
    }


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

        return reduce(source, [:]) { Map map, item ->
                def key = criteria.call(item)
                def value = map.containsKey(key) ? map.get(key)+1 : 1
                map.put(key, value)
                return map
        }
    }

    static public final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> channel) {
        reduce(channel) { min, val -> val<min ? val : min }
    }

    static public final <V> DataflowReadChannel<V> min(final DataflowReadChannel<V> channel, Closure<V> comparator) {

        def _closure
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            _closure = (Closure<V>){ max, item -> comparator.call(item) < comparator.call(max) ? item : max  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            _closure = (Closure<V>){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        reduce(channel, _closure)
    }

    static public final <V> DataflowReadChannel<V>  min(final DataflowQueue<V> channel, Comparator comparator) {
        reduce(channel) { a, b -> comparator.compare(a,b)<0 ? a : b }
    }

    static public final <V> DataflowReadChannel<V> max(final DataflowQueue channel) {
        reduce(channel) { max, val -> val>max ? val : max }
    }

    static public final <V> DataflowReadChannel<V> max(final DataflowQueue<V> channel, Closure comparator) {

        def _closure
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            _closure = (Closure<V>){ max, item -> comparator.call(item) > comparator.call(max) ? item : max  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            _closure = (Closure<V>){ a, b ->  comparator.call(a,b)>0 ? a : b  }
        }
        else {
            throw new IllegalArgumentException("Comparator closure can accept at most 2 arguments")
        }

        reduce(channel, _closure)
    }

    static public final <V> DataflowVariable<V> max(final DataflowQueue<V> channel, Comparator<V> comparator) {
        reduce(channel) { a, b -> comparator.compare(a,b)>0 ? a : b }
    }

    static public final <V> DataflowReadChannel<V> sum(final DataflowQueue<V> channel) {
        reduce(channel, 0) { sum, val -> sum += val }
    }

    static public final <V> DataflowReadChannel<V> sum(final DataflowQueue<V> channel, Closure<V> closure) {
        reduce(channel, 0) { sum, val -> sum += closure.call(val) }
    }


    /**
     * Sorts all collection members into groups determined by the supplied mapping closure
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final DataflowReadChannel<Map> groupBy(final DataflowQueue channel, final Closure criteria ) {

        return reduce(channel, [:]) { map, item ->
            def key = criteria.call(item)
            def list = map.get(key)
            list = list ? list << item : [item]
            map.put(key, list)
            return map
        }
    }

    static public final List<DataflowQueue> separate( final DataflowQueue channel, int n, Closure mapper ) {

        def outputs = []
        n.times { outputs << new DataflowQueue() }
        Dataflow.operator([channel], outputs, new SeparationClosure(mapper))
        outputs
    }

    static public final DataflowQueue spread( final DataflowQueue channel, DataflowReadChannel other ) {

        DataflowQueue result = new DataflowQueue()

        def source
        switch(other) {
            case DataflowQueue:  source = ((DataflowQueue) other).toList(); break
            case DataflowVariable: source = other; break
            default: throw new IllegalArgumentException()
        }

        Dataflow.operator( [channel, source], [result] ) { a, b ->
            [ [a], (b as List) ]
                    .combinations()
                    .each{ Collection it -> bindOutput(it.flatten())  }
        }

        return result
    }


    static public final DataflowQueue spread( final DataflowQueue channel, Collection other ) {
        assert other != null

        return spread(channel, Channel.just(other))
    }


    static public final DataflowReadChannel flatten( final DataflowReadChannel source )  {

        DataflowQueue target = new DataflowQueue()

        def listeners = []

        if( source instanceof DataflowExpression ) {
            listeners << new DataflowEventAdapter() {
                @Override
                public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                    processor.bindOutput( Channel.STOP )
                    processor.terminate()
                }
            }
        }


        Dataflow.operator(inputs: [source], outputs: [target], listeners: listeners) {  it ->
            switch( it ) {
                case Iterable:
                    it.each { value -> bindOutput(value) }
                    break

                case Object[]:
                    it.each { value -> bindOutput(value) }
                    break

                case Map:
                    ((Map)it).each { entry -> bindOutput(entry) }
                    break

                default:
                    bindOutput(it)
            }
        }
        return target
    }

    static public final DataflowQueue window( final DataflowQueue channel, Object closingCriteria ) {

        def closure = new BooleanReturningMethodInvoker("isCase");
        return windowImpl(channel, null, { Object it -> closure.invoke(closingCriteria, it) })

    }

    static public final DataflowQueue window( final DataflowQueue channel, Object startingCriteria, Object closingCriteria  ) {
        assert startingCriteria != null
        assert closingCriteria != null

        def c1 = new BooleanReturningMethodInvoker("isCase");
        def c2 = new BooleanReturningMethodInvoker("isCase");

        return windowImpl(channel, {Object it -> c1.invoke(startingCriteria, it)}, {Object it -> c2.invoke(closingCriteria, it)})

    }

    static public final DataflowQueue window( final DataflowQueue channel, Map params ) {

        if( params.count ) {
            windowWithSizeConstraint( channel, (int)params.count, (int)params?.skip ?: 0 )
        }
        else {
            throw new IllegalArgumentException()
        }
    }



    static private windowWithSizeConstraint( final DataflowQueue channel, int size, int skip = 0 ) {
        assert size>0

        def skipCount = 0
        def itemCount = 0

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

        return windowImpl(channel, skip>0 ? startRule : null, closeRule )
    }


    static private final DataflowQueue windowImpl( final DataflowQueue channel, Closure startingCriteria, Closure closeCriteria ) {
        assert closeCriteria

        // the list holding temporary collected elements
        def buffer = []

        // the result queue
        DataflowQueue result = new DataflowQueue();

        // open frame flag
        boolean isOpen = startingCriteria == null

        // the operator collecting the elements
        Dataflow.operator(channel, result) {
            if( isOpen ) {
                buffer << it
            }
            else if( startingCriteria.call(it) ) {
                isOpen = true
                buffer << it
            }

            if( closeCriteria.call(it) ) {
                bindOutput(buffer);
                buffer = []
                // when a *startingCriteria* is defined, close the open frame flag
                isOpen = (startingCriteria == null)
            }

        }

        return result
    }

    /**
     * Similar to https://github.com/Netflix/RxJava/wiki/Combining-Observables#merge
     *
     * @param source
     * @param target
     * @return
     */
    static final DataflowQueue mix( DataflowQueue source, DataflowQueue... target ) {
        assert target.size()>0

        def result = new DataflowQueue()
        def count = new AtomicInteger( target.size()+1 )
        def handlers = [
                onNext: { result << it },
                onComplete: { if(count.decrementAndGet()==0) { result << Channel.STOP } }
        ]

        source.subscribe(handlers)
        target.each{ it.subscribe(handlers) }

        return result
    }

    /**
     * Phase channels
     *
     * @param source
     * @param target
     * @param mapper
     * @return
     */
    static final DataflowQueue phase( DataflowQueue source, DataflowQueue target, Closure mapper ) {

        def result = new DataflowQueue()
        def state = [:]

        source.subscribe( phaseHandlers(state, 2, source, result, mapper) )
        target.subscribe( phaseHandlers(state, 2, target, result, mapper) )

        return result
    }

    /**
     * Phase channels
     *
     * @param source
     * @param target
     * @return
     */
    static final DataflowQueue phase(DataflowQueue source, DataflowQueue target) {
        phase(source, target, { phaseDefaultMapper(it) } )
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
    static private phaseDefaultMapper( def obj )  {

        switch( obj ) {
            case Map:
                def values = ((Map)obj).keySet()
                return values.size() ? values.getAt(0) : null

            case Map.Entry:
                def entry = (Map.Entry) obj
                return entry.key

            case Collection:
                def values = (Collection)obj
                return values.size() ? values.getAt(0) : null

            case Object[]:
                def values = (Object[])obj
                return values.size() ? values[0] : null

            default:
                return obj
        }

    }

    /**
     * Returns the methods {@code OnNext} and {@code onComplete} which will implement the phase logic
     *
     * @param buffer The shared state buffering the channel received values
     * @param count The overall number of channel
     * @param current The current channel
     * @param channel The channel over which the results are sent
     * @param mapper A closure mapping a value to its key
     * @return A map with {@code OnNext} and {@code onComplete} methods entries
     */
    static private final Map phaseHandlers( Map<Object,Map<DataflowReadChannel,List>> buffer, int count, DataflowReadChannel current, DataflowWriteChannel channel, Closure mapper ) {

        [
                onNext: {
                    synchronized (buffer) {  // phaseImpl is NOT thread safe, synchronize it !
                        def entries = phaseImpl(buffer, count, current, it, mapper)
                        if( entries ) channel.bind(entries)
                    }},

                onComplete: { channel.bind(Channel.STOP) }

        ]

    }

    /**
     * Implements the phase operator logic. Basically buffers the values received on each channel by their key .
     *
     * When a value with the same key has arrived on each channel, they are removed from the buffer and returned as list
     *
     *
     * @param buffer The shared state buffer
     * @param count The overall number of channels
     * @param current The current channel
     * @param item The value just arrived
     * @param mapper The mapping closure retrieving a key by the item just arrived over the current channel
     * @return The list a values having a common key for each channel or {@code null} when some values are missing
     *
     */
    static private final List phaseImpl( Map<Object,Map<DataflowReadChannel,List>> buffer, int count, DataflowReadChannel current, def item, Closure mapper ) {

        // get the index key for this object
        def key = mapper.call(item)

        // given a key we expect to receive on object with the same key on each channel
        def channels = buffer.get(key)
        if( channels==null ) {
            channels = [:]
            buffer[key] = channels
        }

        def entries = channels[current]
        if( entries==null ) {
            entries = []
            channels[current] = entries
        }

        // add the received item to the list
        entries << item

        // now check if it has received a element matching for each channel
        if( channels.size() != count )  {
            return null
        }

        def result = []

        Iterator<Map.Entry<DataflowReadChannel,List>> itr = channels.iterator()
        while( itr.hasNext() ) {
            def entry = itr.next()

            def list = entry.getValue()
            result << list[0]
            list.remove(0)
            if( list.size() == 0 ) {
                itr.remove()
            }
        }

        return result
    }

    /**
     * Makes the output of the source channel to be an input for the specified channels
     *
     * @param source The source dataflow object
     * @param target One or more writable to which source is copied
     */
    public static <T> void split( DataflowReadChannel<T> source, DataflowWriteChannel<T>... target ) {
        assert source != null
        assert target

        source.split( target as List )
    }

}
