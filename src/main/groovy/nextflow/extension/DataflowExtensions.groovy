package nextflow.extension
import static java.util.Arrays.asList

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.SeparationClosure
import nextflow.Channel
import nextflow.util.Duration
import org.codehaus.groovy.runtime.NullObject
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class DataflowExtensions {

    /**
     * Subscribe *onNext* event
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> subscribe(final DataflowQueue channel, final Closure<V> closure) {
        subscribe( channel, [onNext: closure] )
    }


    /**
     * Subscribe *onNext*, *onError* and *onComplete*
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowReadChannel<V> subscribe(final DataflowQueue channel, final Map<String,Closure> events ) {

        if( !events.onNext && !events.onComplete && !events.onError ) {
            throw new IllegalArgumentException("You must specify at least an event between: onNext, onComplete, onError")
        }

        def listener = new DataflowEventAdapter() {

            int index = 0

            @Override
            public void afterRun(final DataflowProcessor processor, final List<Object> messages) {

                def closure = events.onNext
                if( closure ) {
                    if( closure.getMaximumNumberOfParameters() == 1 ) closure.call( messages.get(0) )
                    if( closure.getMaximumNumberOfParameters() >= 2 ) closure.call( messages.get(0), this.index++ )
                }
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {

                def closure = events.onError
                if( closure ) {
                    DataflowExtensions.log.debug('Dataflow subscribe() exception', e)

                    def result = null
                    if( closure.getMaximumNumberOfParameters() == 1 ) result = closure.call( e )
                    if( closure.getMaximumNumberOfParameters() >= 2 ) result = closure.call( e, this.index )

                    if( result instanceof Boolean ) {
                        return result;
                    }
                }
                else {
                    DataflowExtensions.log.error('Dataflow subscribe() exception', e)
                }

                return true
            }

            @Override
            public void afterStop(final DataflowProcessor processor) {

                def closure = events.onComplete
                if( closure ) {
                    closure.call()
                }

            }
        }

        final DataflowQueue<V> result = new DataflowQueue<V>();
        final Map<String, Object> parameters = new HashMap<String, Object>();
        parameters.put("inputs", [channel]);
        parameters.put("outputs", [result]);
        parameters.put('listeners', [listener])

        Dataflow.retrieveCurrentDFPGroup().operator(parameters, {it});

        // forwards all
        return channel

    }


    /**
     * Just a synonym for {@code DataflowQueue#chainWith} method
     *
     * @param channel
     * @param closure
     * @return
     */
    static public final <V> DataflowQueue<V> map(final DataflowQueue channel, final Closure<V> closure) {
        return channel.chain(closure)
    }

    static public final <V> DataflowQueue<V> mapMany(final DataflowQueue channel, final Closure<V> closure) {

        def result = new DataflowQueue<V>()

        Dataflow.operator(channel, result) {

            def value = closure.call(it);
            if( value instanceof Collection ) {
                value.each{ entry -> bindOutput(entry) }
            }
            else {
                bindOutput(value)
            }
        }

        return result
    }

    static public final <V> DataflowVariable<V> reduce(final DataflowQueue<V> channel, final Closure<V> closure) {

        // the dataflow variable to return the final aggregation value
        def result = new DataflowVariable()

        // the *accumulator* value
        def count = 0
        def accum = null

        // intercepts operator events
        def listener = new DataflowEventAdapter() {
            /*
             * call the passed closure each time
             */
            public void afterRun(final DataflowProcessor processor, final List<Object> messages) {
                def item = messages.get(0)
                def value
                if( count++ == 0 ) {
                    value = item
                }
                else {
                    value = closure.call(accum, item)
                }

                if( value == Channel.STOP ) {
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
     *
     *
     * @param channel
     * @param seed
     * @param closure
     * @return
     */
    static public final <V> DataflowVariable<V> reduce(final DataflowQueue channel, V seed, final Closure<V> closure) {

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
                final value = closure.call(accum, item)
                if( value == Channel.STOP ) {
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
    static public final <V> DataflowQueue<V> grep(final DataflowQueue channel, final Object criteria) {
        def discriminator = new BooleanReturningMethodInvoker("isCase");
        return channel.filter { Object it -> discriminator.invoke(criteria, it) }
    }

    public static <V> DataflowReadChannel<V> chain(final DataflowQueue<V> channel, final Closure<V> closure) {
        final DataflowQueue<V> result = new DataflowQueue<V>();
        Dataflow.operator(channel, result, new ChainWithClosure<V>(closure));
        return result;
    }

    public static <V> DataflowReadChannel<V> chain(final DataflowQueue<V> channel, final Map<String, Object> params, final Closure<V> closure) {
        final DataflowQueue<V> result = new DataflowQueue<V>();
        final Map<String, Object> parameters = new HashMap<String, Object>(params);
        parameters.put("inputs", asList(channel));
        parameters.put("outputs", asList(result));

        Dataflow.operator(parameters, new ChainWithClosure<V>(closure));
        return result;
    }

    /**
     * Modifies this collection to remove all duplicated items, using the default comparator.
     *
     * assert [1,3] == [1,3,3].unique()
     *
     * @param channel
     * @return
     */
    static public final <V> DataflowQueue<V> unique(final DataflowQueue<V> channel) {
        unique(channel) { it }
    }

    /**
     * A convenience method for making a collection unique using a Closure to determine duplicate (equal) items. If the closure takes a single parameter, the argument passed will be each element, and the closure should return a value used for comparison (either using Comparable#compareTo or Object#equals). If the closure takes two parameters, two items from the collection will be passed as arguments, and the closure should return an int value (with 0 indicating the items are not unique).
     * assert [1,4] == [1,3,4,5].unique { it % 2 }
     * assert [2,3,4] == [2,3,3,4].unique { a, b -> a <=> b }
     *
     * @param channel
     * @param comparator
     * @return
     */
    static public final <V> DataflowQueue<V> unique(final DataflowQueue channel, Closure comparator ) {

        def history = [:]

        // when the operator stop clear the history map
        def events = new DataflowEventAdapter() {
            public void afterStop(final DataflowProcessor processor) {
                history.clear()
            }
        }

        def filter = {
            def key = comparator.call(it)
            if( history.containsKey(key) ) {
                return NullObject.getNullObject()
            }
            else {
                history.put(key,true)
                return it
            }
        }  as Closure<V>

        // filter removing all duplicates
        return (DataflowQueue<V>)channel.chain(listeners: [events], filter)

    }

    /**
     * Suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     *
     * @return
     */
    static public final <V> DataflowQueue<V> distinct( final DataflowQueue<V> channel ) {
        distinct(channel) {it}
    }

    /**
     * suppress duplicate consecutive items emitted by the source Observable
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#suppress-duplicate-consecutive-items-emitted-by-the-source-observable
     *
     * @return
     */
    static public final <V> DataflowQueue<V> distinct( final DataflowQueue channel, Closure comparator ) {

        V previous = null

        return channel.chain {

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
     * Emit only the first item emitted by an Observable, or the first item that meets some condition
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#first
     *
     * @param channel
     * @return
     */
    static public final <V> DataflowVariable<V> first( final DataflowQueue<V> channel ) {

        def result = new DataflowVariable<V>()
        channel.whenBound { result.bind(it) }
        return result
    }

    /**
     *
     * Emit only the first item emitted by an Observable, or the first item that meets some condition
     *
     * See https://github.com/Netflix/RxJava/wiki/Filtering-Observables#first
     *
     * @param channel
     * @return
     */
    static public final <V> DataflowVariable<V> first( final DataflowQueue channel, Closure criteria ) {

        def result = new DataflowVariable()
        Dataflow.operator([channel],[]) {

            if( criteria.call(it)) {
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
     * @param channel
     * @return
     */
    static public final <V> DataflowQueue<V> take( final DataflowQueue<V> channel, int n ) {

        def count = 0
        def result = new DataflowQueue<V>()
        Dataflow.operator([channel],[]) {

            if( count++ < n ) {
                result.bind(it)
                return
            }
            result.bind(Channel.STOP)
            terminate()
        }

        return result
    }


    static public final <V> DataflowVariable<V> last( final DataflowQueue channel  ) {

        def result = new DataflowVariable()
        def last = null
        channel.subscribe( onNext: { last = it }, onComplete: {  result.bind(last) } )
        return result

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
    static public final <V> DataflowQueue<V> parallel( final DataflowQueue<V> channel, Closure closure ) {

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
    static public final <V> DataflowReadChannel<V> doFinally( final DataflowQueue channel, Closure<V> closure ) {
        channel.subscribe( onComplete: closure )
    }

    /**
     * Convert a {@code DataflowQueue} alias *channel* to a Java {@code List}
     *
     * @param channel The channel to be converted
     * @return A list holding all the items send over the channel
     */
    static public final <V> DataflowReadChannel<V> toList(final DataflowQueue channel) {
        return reduce(channel, []) { list, item -> list << item }
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param channel
     * @param value
     * @return
     */
    static public final DataflowVariable<Number> count(final DataflowQueue channel ) {
        reduce(channel, 0) { current, item -> current+1 }
    }

    /**
     * Counts the number of occurrences of the given value inside this collection.
     *
     * @param channel
     * @param value
     * @return
     */
    static public final DataflowVariable<Number> count(final DataflowQueue channel, final Object value ) {
        reduce(channel, 0) { current, item -> item==value ? current+1 : current }
    }

    /**
     * Counts the number of occurrences which satisfy the given closure from inside this collection
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final DataflowVariable<Number> count(final DataflowQueue channel, final Closure<Boolean> criteria ) {
        reduce(channel, 0) { current, item -> criteria.call(item) ? current+1 : current }
    }

    /**
     * Sorts all collection members into groups determined by the supplied mapping closure and counts the group size
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final DataflowVariable<Map> countBy(final DataflowQueue channel, final Closure criteria ) {

        return reduce(channel, [:]) { Map map, item ->
                def key = criteria.call(item)
                def value = map.containsKey(key) ? map.get(key)+1 : 1
                map.put(key, value)
                return map
        }
    }

    static public final <V> DataflowVariable<V> min(final DataflowQueue<V> channel) {
        reduce(channel) { min, val -> val<min ? val : min }
    }

    static public final <V> DataflowVariable<V> min(final DataflowQueue<V> channel, Closure<V> comparator) {

        def _closure
        if( comparator.getMaximumNumberOfParameters() == 1 ) {
            _closure = (Closure<V>){ max, item -> comparator.call(item) < comparator.call(max) ? item : max  }
        }
        else if( comparator.getMaximumNumberOfParameters() == 2 ) {
            _closure = (Closure<V>){ a, b ->  comparator.call(a,b) < 0 ? a : b  }
        }

        reduce(channel, _closure)
    }

    static public final <V> DataflowVariable<V>  min(final DataflowQueue channel, Comparator comparator) {
        reduce(channel) { a, b -> comparator.compare(a,b)<0 ? a : b }
    }

    static public final <V> DataflowVariable<V> max(final DataflowQueue channel) {
        reduce(channel) { max, val -> val>max ? val : max }
    }

    static public final <V> DataflowVariable<V> max(final DataflowQueue<V> channel, Closure comparator) {

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

    static public final <V> DataflowVariable<V> sum(final DataflowQueue<V> channel) {
        reduce(channel, 0) { sum, val -> sum += val }
    }

    static public final <V> DataflowVariable<V> sum(final DataflowQueue<V> channel, Closure<V> closure) {
        reduce(channel, 0) { sum, val -> sum += closure.call(val) }
    }


    /**
     * Sorts all collection members into groups determined by the supplied mapping closure
     *
     * @param channel
     * @param criteria
     * @return
     */
    static public final DataflowVariable<Map> groupBy(final DataflowQueue channel, final Closure criteria ) {

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


    static public final DataflowQueue flatten( final DataflowQueue channel )  {

        DataflowQueue result = new DataflowQueue()

        Dataflow.operator(channel,result) {  it ->
            if( it instanceof Collection ) { it.each { value -> bindOutput(value) } }
            else { bindOutput(it) }
        }
        return result
    }

    static public final DataflowQueue buffer( final DataflowQueue channel, Object closingCriteria ) {

    }

    static public final DataflowQueue buffer( final DataflowQueue channel, long count ) {

    }

    static public final DataflowQueue buffer( final DataflowQueue channel, long count, long skip ) {

    }


    static public final DataflowQueue buffer( final DataflowQueue channel, Duration timespan ) {

    }

    static public final DataflowQueue buffer( final DataflowQueue channel, Duration timespan, long count ) {

    }


    static public final DataflowQueue buffer( final DataflowQueue channel, Object startingCriteria, Object closingCriteria  ) {

    }
}
