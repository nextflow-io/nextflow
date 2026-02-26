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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.extension.op.ContextRunPerThread
import nextflow.extension.op.Op
import nextflow.extension.op.OpContext
import nextflow.extension.op.OpDatum
import nextflow.prov.OperatorRun
import nextflow.util.ArrayBag
import nextflow.util.CacheHelper
import nextflow.util.CheckHelper
/**
 * Implements {@link OperatorImpl#groupTuple} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GroupTupleOp {

    static final private Map GROUP_TUPLE_PARAMS = [ by: [Integer, List], sort: [Boolean, 'true','natural','deep','hash',Closure,Comparator], size: Integer, remainder: Boolean ]

    static final private List<Integer> GROUP_DEFAULT_INDEX = [0]

    /**
     * Comparator used to sort tuple entries (when required)
     */
    private Comparator comparator

    private Integer size

    private List indices

    private DataflowWriteChannel target

    private DataflowReadChannel channel

    private boolean remainder

    private Map<List,List> groups = [:]

    private sort

    private OpContext context = new ContextRunPerThread()

    GroupTupleOp(Map params, DataflowReadChannel source) {
        CheckHelper.checkParams('groupTuple', params, GROUP_TUPLE_PARAMS)

        channel = source
        indices = getGroupTupleIndices(params)
        size = params?.size as Integer ?: 0
        remainder = params?.remainder ?: false
        sort = params?.sort

        defineComparator()
    }

    GroupTupleOp withTarget(DataflowWriteChannel target) {
        this.target = target
        return this
    }

    static private List<Integer> getGroupTupleIndices( Map params ) {

        if( params?.by == null )
            return GROUP_DEFAULT_INDEX

        if( params.by instanceof List )
            return params.by as List<Integer>

        if( params.by instanceof Integer || params.by.toString().isInteger() )
            return [params.by as Integer]

        throw new IllegalArgumentException("Not a valid `by` index for `groupTuple` operator: '${params.by}' -- It must be an integer value or a list of integers")
    }

    /*
     * Collects received values grouping by key
     */
    private void collectTuple(List tuple) {

        final key = tuple[indices]                      // the actual grouping key
        final len = tuple.size()

        final List items = groups.getOrCreate(key) {    // get the group for the specified key
            def result = new ArrayList(len)             // create if does not exist
            for( int i=0; i<len; i++ )
                result[i] = (i in indices ? tuple[i] : new ArrayBag())
            return result
        }

        final run = context.getOperatorRun()
        int count=-1
        for( int i=0; i<len; i++ ) {                    // append the values in the tuple
            if( i !in indices ) {
                def list = (items[i] as List)
                if( list==null ) {
                    list = new ArrayBag()
                    items.add(i, list)
                }
                // wrap the acquired value in OpDatum object to track the input provenance
                list.add( OpDatum.of(tuple[i], run) )
                count = list.size()
            }
        }

        final sz = size ?: sizeBy(key)
        if( sz>0 && sz==count ) {
            bindTuple(items, sz)
            groups.remove(key)
        }
    }

    /*
     * finalize the grouping binding the remaining values
     */
    private void finalise(DataflowProcessor dp) {
        groups.each { keys, items -> bindTuple(items, size ?: sizeBy(keys)) }
        Op.bind(dp, target, Channel.STOP)
    }

    /*
     * bind collected items to the target channel
     */
    private void bindTuple( List items, int sz ) {
        final tuple = new ArrayList(items)
        if( !remainder && sz>0 ) {
            // verify exist it contains 'size' elements
            def list = (List) items.find { it instanceof List }
            if( list.size() != sz ) {
                return
            }
        }
        // unwrap all "OpData" object and restore original values
        final run = unwrapValues(tuple)
        // sort the tuple content when a comparator is defined
        if( comparator ) {
            sortInnerLists(tuple, comparator)
        }
        // finally bind the resulting tuple
        Op.bind(run, target, tuple)
    }

    static protected OperatorRun unwrapValues(List tuple) {
        final inputs = new ArrayList<Integer>()

        for( Object it : tuple ) {
            if( it instanceof ArrayBag ) {
                final bag = it
                for( int i=0; i<bag.size(); i++ ) {
                    bag[i] = OpDatum.unwrap(bag[i], inputs)
                }
            }
        }

        return new OperatorRun(new LinkedHashSet<Integer>(inputs))
    }

    /**
     * Define the comparator to be used depending the #sort property
     */
    private void defineComparator( ) {
        /*
         * comparator logic used to sort tuple elements
         */
        switch(sort) {
            case null:
                break

            case true:
            case 'true':
            case 'natural':
                comparator = { o1, o2 -> o1<=>o2 } as Comparator<Comparable>
                break;

            case 'hash':
                comparator = { o1, o2 ->
                    final h1 = CacheHelper.hasher(o1).hash()
                    final h2 = CacheHelper.hasher(o2).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator
                break

            case 'deep':
                comparator = { o1, o2 ->
                    final h1 = CacheHelper.hasher(o1, CacheHelper.HashMode.DEEP).hash()
                    final h2 = CacheHelper.hasher(o2, CacheHelper.HashMode.DEEP).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator
                break

            case Comparator:
                comparator = sort as Comparator
                break

            case Closure:
                final closure = (Closure)sort
                if( closure.getMaximumNumberOfParameters()==2 ) {
                    comparator = sort as Comparator
                }
                else if( closure.getMaximumNumberOfParameters()==1 ) {
                    comparator = { o1, o2 ->
                        final v1 = closure.call(o1) as Comparable
                        final v2 = closure.call(o2) as Comparable
                        return v1 <=> v2
                    } as Comparator
                }
                else
                    throw new IllegalArgumentException("Invalid groupTuple option - The closure should have 1 or 2 arguments")
                break

            default:
                throw new IllegalArgumentException("Not a valid sort argument: ${sort}")
        }
    }

    /**
     * Main method to invoke the operator
     *
     * @return The resulting channel
     */
    DataflowWriteChannel apply() {

        if( target == null )
            target = CH.create()

        /*
         * apply the logic to the source channel
         */
        new SubscribeOp()
            .withInput(channel)
            .withContext(context)
            .withOnNext(this.&collectTuple)
            .withOnComplete(this.&finalise)
            .apply()

        /*
         * return the target channel
         */
        return target
    }

    private static sortInnerLists( List tuple, Comparator c ) {
        for( int i=0; i<tuple.size(); i++ ) {
            final entry = tuple[i]
            if( entry !instanceof List ) continue
            Collections.sort(entry as List, c)
        }
    }

    static protected int sizeBy(List target)  {
        if( target.size()==1 && target[0] instanceof GroupKey ) {
            final group = (GroupKey)target[0]
            final size = group.getGroupSize()
            log.debug "GroupTuple dynamic size: key=${group} size=$size"
            return size
        }
        else
            return 0
    }

}
