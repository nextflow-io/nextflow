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

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.NF
import nextflow.exception.AbortOperationException
import nextflow.util.CheckHelper
/**
 * Implements {@link OperatorImpl#join} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class JoinOp {

    static final private Map JOIN_PARAMS = [remainder: Boolean, by: [List,Integer], failOnMismatch: Boolean, failOnDuplicate: Boolean]

    private DataflowReadChannel source

    private DataflowReadChannel target

    private List<Integer> pivot

    private boolean remainder

    private Boolean[] singleton = [null,null] as Boolean[]

    private boolean failOnMismatch

    private boolean failOnDuplicate

    private volatile boolean failed

    private Set uniqueKeys = new LinkedHashSet()

    JoinOp( DataflowReadChannel source, DataflowReadChannel target, Map params = null ) {
        CheckHelper.checkParams('join', params, JOIN_PARAMS)
        this.source = source
        this.target = target
        this.pivot = parsePivot(params?.by)
        this.remainder = params?.remainder ? params.remainder as boolean : false
        this.failOnMismatch = params?.failOnMismatch || (!this.remainder && NF.isStrictMode())
        this.failOnDuplicate = params?.failOnDuplicate || NF.isStrictMode()
        // sanity check
        if( remainder && failOnMismatch )
            throw new AbortOperationException("Conflicting join operator options `remainder` and `failOnMismatch`")
    }

    private List<Integer> parsePivot(value) {
        if( value == null )
            return [0]
        if( value instanceof List )
            return value as List<Integer>
        else
            return [value as int]
    }

    DataflowWriteChannel apply() {

        // the resulting channel
        final result = CH.create()
        // the following buffer maintains the state of collected items as a map of maps.
        // The first map associates the joining key with the collected values
        // The inner map associates the channel index with the actual values received on that channel
        final Map<Object,Map<Integer,List>> state = [:]

        final count = 2
        final stopCount = new AtomicInteger(count)

        DataflowHelper.subscribeImpl( source, handler(state, count, 0, result, stopCount, remainder) )
        DataflowHelper.subscribeImpl( target, handler(state, count, 1, result, stopCount, remainder) )
        return result
    }

    /**
     * Returns the methods {@code OnNext} and {@code onComplete} which will implement the join logic
     *
     * @param buffer The shared state buffering the channel received values
     * @param count The overall number of channel
     * @param current The current channel
     * @param target The channel over which the results are sent
     * @param mapper A closure mapping a value to its key
     * @return A map with {@code OnNext} and {@code onComplete} methods entries
     */
    private final Map<String,Closure> handler( Map<Object,Map<Integer,List>> buffer, int size, int index, DataflowWriteChannel target, AtomicInteger stopCount, boolean remainder ) {

        final Map<String,Closure> result = new HashMap<>(2)

        result.onNext = {
            synchronized (this) {
                if(!failed) try {
                    def entries = join0(buffer, size, index, it)
                    if( entries ) {
                        target.bind( entries.size()==1 ? entries[0] : entries )
                    }
                }
                catch (Exception e) {
                    failed = true
                    target << Channel.STOP
                    throw e
                }
            }}

        result.onComplete = {
            if( stopCount.decrementAndGet()==0 && !failed ) {
                try {
                    if( remainder || failOnDuplicate )
                        checkRemainder(buffer,size,target)
                    if( failOnMismatch )
                        checkForMismatch(buffer)
                }
                finally {
                    target << Channel.STOP
                }
            }}
        
        return result
    }


    /**
     * Implements the join operator logic. Basically buffers the values received on each channel by their key .
     *
     * When a value with the same key has arrived on each channel, they are removed from the buffer and returned as list
     *
     *
     * @param buffer The shared state buffer
     * @param size The overall number of channels
     * @param current The current channel
     * @param data The value just arrived
     * @param mapper The mapping closure retrieving a key by the item just arrived over the current channel
     * @return The list a values having a common key for each channel or {@code null} when some values are missing
     *
     */
    @PackageScope
    final List join0( Map<Object,Map<Integer,List>> buffer, int size, int index, def data) {

        // The 'buffer' structure track the values emitted by the channel, it is arranged in the following manner:
        //
        //  Map< key, Map< channel index, List[ values ] >  >
        //
        // In the main map there's an entry for each 'key' for which a match is required,
        // to which is associated another map which associate the channel (index) on which the item
        // has been emitted and all the values received (for that channel) not yet emitted.
        // (this is required to do not lost an item that is emitted more than one time on the same channel
        //  before a match for it is found on another channel)

        // get the index key for this object
        final item0 = DataflowHelper.makeKey(pivot, data)

        // check for unique keys
        checkForDuplicate(item0.keys, item0.values, index, false)

        // given a key we expect to receive on object with the same key on each channel
        Map<Integer,List> channels = buffer.get(item0.keys)
        if( channels==null ) {
            channels = new TreeMap<Integer, List>()
            buffer[item0.keys] = channels
        }

        if( !channels.containsKey(index) ) {
            channels[index] = []
        }

        def entries = channels[index]

        // add the received item to the list
        // when it is used in the gather op add always as the first item
        entries << item0.values
        setSingleton(index, item0.values.size()==0)

        // now check if it has received an element matching for each channel
        if( channels.size() != size )  {
            return null
        }

        def result = []
        // add the key
        addToList(result, item0.keys)

        final itr = channels.iterator()
        while( itr.hasNext() ) {
            def entry = (Map.Entry<Integer,List>)itr.next()

            def list = entry.getValue()
            addToList(result, list[0])

            list.remove(0)
            if( list.size() == 0 ) {
                itr.remove()
            }
        }

        // track unique keys
        checkForDuplicate(item0.keys, item0.values, index, true)
        return result
    }

    private final void checkRemainder(Map<Object,Map<Integer,List>> buffers, int count, DataflowWriteChannel target ) {
       log.trace "Operator `join` remainder buffer: ${-> buffers}"

        for( Object key : buffers.keySet() ) {
            Map<Integer,List> entry = buffers.get(key)

            while( true ) {

                boolean fill=false
                def result = new ArrayList(count+1)
                addToList(result, key)

                for( int i=0; i<count; i++ ) {
                    List values = entry[i]
                    if( values ) {
                        // track unique keys
                        checkForDuplicate(key, values[0], i,false)

                        addToList(result, values[0])
                        fill |= true
                        values.remove(0)
                    }
                    else if( !singleton(i) ) {
                        addToList(result, null)
                    }
                }

                if( fill ) {
                    final value = singleton() ? result[0] : result
                    // bind value to target channel
                    if( remainder ) target.bind(value)
                }
                else
                    break
            }

        }
    }

    protected void checkForDuplicate( key, value, int dir, boolean add ) {
        if( failOnDuplicate && ( (add && !uniqueKeys.add(key)) || (!add && uniqueKeys.contains(key)) ) ) {
            final msg = "Detected join operation duplicate emission on ${dir==0 ? 'left' : 'right'} channel -- offending element: key=${csv0(key,',')}; value=${csv0(value,',')}"
            throw new AbortOperationException(msg)
        }
    }

    protected void checkForMismatch( Map<Object,Map<Integer,List>> buffers ) {
        final result = new HashMap<Object,List>()
        for( Object key : buffers.keySet() ) {
            Map<Integer,List> el = buffers.get(key)
            final remainder0 = el.entrySet()
            if( !remainder0 )
                continue

            result[key] = []
            for( Map.Entry entry : remainder0 ) {
                result[key].add(csv0(entry.value,','))
            }
        }

        if( result ) {
            final list = result.take(10).collect { k,v -> "key=${csv0(k,';')} values=${csv0(v,';')}" }
            final str = result.size() == 1 ? list[0] : list.collect { "\n- $it" }.join(' ')
            def message = "Join mismatch for the following entries: $str"
            if( list.size()!=result.size() )
                message += '\n(more omitted)'
            throw new AbortOperationException(message)
        }
    }


    private String csv0(value, String sep) {
        value instanceof List ? value.join(sep) : value.toString()
    }

    private boolean singleton(int i=-1) {
        if(i==-1)
            return singleton[0]!=false && singleton[1]!=false

        def result = singleton[i]
        if( result == null )
            result = singleton[ (i+1)%2 ]

        return result
    }

    private void setSingleton(int index, boolean flag) {
        singleton[index] = singleton[index]==null ? flag : singleton[index] || flag
    }

}
