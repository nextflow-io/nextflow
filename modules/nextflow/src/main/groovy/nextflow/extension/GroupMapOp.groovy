/*
 * Copyright 2023, Seqera Labs
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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.util.ArrayBag
import nextflow.util.CacheHelper
import nextflow.util.CheckHelper
/**
 * Implements {@link OperatorImpl#groupMap} operator logic
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class GroupMapOp extends AbstractGroupOp {

    private List indices

    private Map<List,Map> groups = [:]

    GroupMapOp(Map params, DataflowReadChannel source) {

        super(params, source)

        CheckHelper.checkParams('groupMap', params, PARAM_TYPES + [ by: [String, List] ])

        indices = getIndices(params?.by)
    }

    static private List<String> getIndices( by ) {

        if( by == null )
            throw new IllegalArgumentException("The `by` option is required for `groupMap` operator: '${by}'")

        if( by instanceof List )
            return by

        if( by instanceof String )
            return [by]

        throw new IllegalArgumentException("Not a valid `by` index for `groupMap` operator: '${by}' -- It must be a string or a list of strings")
    }

    @Override
    protected Map getHandlers() {
        [onNext: this.&onNext, onComplete: this.&onComplete]
    }

    /**
     * Collect a received item into its group.
     *
     * @param item
     */
    private void onNext(Map item) {

        // get the grouping key
        final key = indices.collect { k -> item[k] }

        // get the group for the specified key
        // or create it if it does not exist
        final Map group = groups.getOrCreate(key) {
            def result = new HashMap(item.size())
            for( String k : item.keySet() )
                result[k] = (k in indices ? item[k] : new ArrayBag())
            return result
        }

        // append the values in the item
        int count = -1
        for( String k : item.keySet() ) {
            if( k !in indices ) {
                def list = (List)group[k]
                list.add( item[k] )
                count = list.size()
            }
        }

        // emit group if it is complete
        final size = this.size ?: sizeBy(key)
        if( size > 0 && size == count ) {
            bindGroup(group, size)
            groups.remove(key)
        }
    }

    /**
     * Emit the remaining groups when all values have been received.
     */
    private void onComplete(nop) {
        groups.each { key, group -> bindGroup(group, size ?: sizeBy(key)) }
        target.bind(Channel.STOP)
    }

    /**
     * Emit a group.
     *
     * @param group
     * @param size
     */
    private void bindGroup( Map group, int size ) {

        def item = new HashMap(group)

        if( !remainder && size > 0 ) {
            // make sure the group contains 'size' elements
            List list = group.values().find { it instanceof List }
            if( list.size() != size )
                return
        }

        // sort the grouped entries
        if( comparator )
            for( def entry : item.values() )
                if( entry instanceof List )
                    Collections.sort((List)entry, comparator)

        target.bind( item )
    }

}
