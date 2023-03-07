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
 * Abstract class for grouping operators.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
abstract class AbstractGroupOp {

    static protected Map PARAM_TYPES = [
        sort: [Boolean, 'true', 'natural', 'deep', 'hash', Closure, Comparator],
        size: Integer,
        remainder: Boolean
    ]

    protected DataflowReadChannel source

    protected int size

    protected boolean remainder

    protected sort

    protected Comparator comparator

    protected DataflowWriteChannel target

    AbstractGroupOp(Map params, DataflowReadChannel source) {
        this.source = source
        size = params?.size ?: 0
        remainder = params?.remainder ?: false
        sort = params?.sort
        comparator = createComparator()
    }

    /**
     * Get the expected size of a key that contains a GroupKey.
     *
     * @param key
     */
    static protected int sizeBy(List key)  {
        if( key.size()==1 && key[0] instanceof GroupKey ) {
            final groupKey = (GroupKey)key[0]
            final size = groupKey.getGroupSize()
            log.debug "groupMap dynamic size: key=${groupKey} size=${size}"
            return size
        }
        else
            return 0
    }

    /**
     * Create the comparator to be used depending the #sort property
     */
    private Comparator createComparator() {

        /*
         * comparator logic used to sort tuple elements
         */
        switch(sort) {
            case null:
                return null

            case true:
            case 'true':
            case 'natural':
                return { o1,o2 -> o1<=>o2 } as Comparator

            case 'hash':
                return { o1, o2 ->
                    def h1 = CacheHelper.hasher(o1).hash()
                    def h2 = CacheHelper.hasher(o2).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator

            case 'deep':
                return { o1, o2 ->
                    def h1 = CacheHelper.hasher(o1, CacheHelper.HashMode.DEEP).hash()
                    def h2 = CacheHelper.hasher(o2, CacheHelper.HashMode.DEEP).hash()
                    return h1.asLong() <=> h2.asLong()
                } as Comparator

            case Comparator:
                return sort as Comparator

            case Closure:
                final closure = (Closure)sort
                if( closure.getMaximumNumberOfParameters()==2 )
                    return sort as Comparator

                else if( closure.getMaximumNumberOfParameters()==1 )
                    return { o1, o2 ->
                        def v1 = closure.call(o1)
                        def v2 = closure.call(o2)
                        return v1 <=> v2
                    } as Comparator

                else
                    throw new IllegalArgumentException("Invalid groupMap option - The sort closure should have 1 or 2 arguments")

            default:
                throw new IllegalArgumentException("Not a valid sort argument: ${sort}")
        }

    }

    /**
     * Invoke the operator
     */
    DataflowWriteChannel apply() {

        if( target == null )
            target = CH.create()

        // apply the operator to the source channel
        DataflowHelper.subscribeImpl(source, getHandlers())

        // return the target channel
        return target
    }

    abstract protected Map getHandlers()

}
