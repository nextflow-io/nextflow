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

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.util.ArrayBag
import nextflow.util.CacheHelper
import nextflow.util.CheckHelper
/**
 * Implements {@link DataflowExtensions#groupTuple} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GroupTupleOp {

    static private Map GROUP_TUPLE_PARAMS = [ by: [Integer, List], sort: [Boolean, 'true','natural','deep','hash',Closure,Comparator], size: Integer, remainder: Boolean ]

    static private List<Integer> GROUP_DEFAULT_INDEX = [0]


    /**
     * Comparator used to sort tuple entries (when required)
     */
    private Comparator comparator

    private int size

    private List indices

    private DataflowQueue target

    private DataflowReadChannel channel

    private boolean remainder

    private Map<List,List> groups = [:]

    private sort

    GroupTupleOp(Map params, DataflowReadChannel source) {

        CheckHelper.checkParams('groupTuple', params, GROUP_TUPLE_PARAMS)

        channel = source
        target = new DataflowQueue()
        indices = getGroupTupleIndices(params)
        size = params?.size ?: 0
        remainder = params?.remainder ?: false
        sort = params?.sort

        defineComparator()
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
    private void collect(List tuple) {

        final key = tuple[indices]                      // the actual grouping key
        final len = tuple.size()

        final List items = groups.getOrCreate(key) {    // get the group for the specified key
            def result = new ArrayList(len)             // create if does not exists
            for( int i=0; i<len; i++ )
                result[i] = (i in indices ? tuple[i] : new ArrayBag())
            return result
        }

        int count=-1
        for( int i=0; i<len; i++ ) {                    // append the values in the tuple
            if( ! (i in indices) ) {
                def list = (items[i] as List)
                list.add( tuple[i] )
                count=list.size()
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
    private void finalise(nop) {
        groups.each { keys, items -> bindTuple(items, size ?: sizeBy(keys)) }
        target.bind(Channel.STOP)
    }

    /*
     * bind collected items to the target channel
     */
    private void bindTuple( List items, int sz ) {

        def tuple = new ArrayList(items)

        if( !remainder && sz>0 ) {
            // verify exist it contains 'size' elements
            List list = items.find { it instanceof List }
            if( list.size() != sz ) {
                return
            }
        }

        if( comparator ) {
            sortInnerLists(tuple, comparator)
        }

        target.bind( tuple )
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
                comparator = { o1,o2 -> o1<=>o2 } as Comparator
                break;

            case 'hash':
                comparator = { o1, o2 ->
                    def h1 = CacheHelper.hasher(o1).hash()
                    def h2 = CacheHelper.hasher(o2).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator
                break

            case 'deep':
                comparator = { o1, o2 ->
                    def h1 = CacheHelper.hasher(o1, CacheHelper.HashMode.DEEP).hash()
                    def h2 = CacheHelper.hasher(o2, CacheHelper.HashMode.DEEP).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator
                break

            case Comparator:
                comparator = sort as Comparator
                break

            case Closure:
                comparator = { o1, o2 ->
                    def closure = (Closure)sort
                    def v1 = closure.call(o1)
                    def v2 = closure.call(o2)
                    return v1 <=> v2
                } as Comparator
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
    DataflowQueue apply() {

        /*
         * apply the logic the the source channel
         */
        DataflowHelper.subscribeImpl(channel, [onNext: this.&collect, onComplete: this.&finalise])

        /*
         * return the target channel
         */
        return target
    }

    private static sortInnerLists( List tuple, Comparator c ) {

        for( int i=0; i<tuple.size(); i++ ) {
            def entry = tuple[i]
            if( !(entry instanceof List) ) continue
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
