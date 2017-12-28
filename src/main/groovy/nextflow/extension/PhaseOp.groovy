/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.extension
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.util.CheckHelper

/**
 * Implements {@link DataflowExtensions#phase} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PhaseOp {

    static final private Map PHASE_PARAMS = [remainder: Boolean]

    private DataflowReadChannel source

    private Map opts

    private DataflowReadChannel target

    private Closure mapper = DataflowExtensions.DEFAULT_MAPPING_CLOSURE

    PhaseOp( DataflowReadChannel source, DataflowReadChannel target ) {
        this.source = source
        this.target = target
        this.opts = [:]
    }

    PhaseOp setOpts( Map opts ) {
        CheckHelper.checkParams('phase', opts, PHASE_PARAMS)
        this.opts = opts
        return this
    }

    PhaseOp setMapper( Closure mapper ) {
        this.mapper = mapper ?: DataflowExtensions.DEFAULT_MAPPING_CLOSURE
        return this
    }

    DataflowQueue apply() {

        def result = new DataflowQueue()
        def state = [:]

        final count = 2
        final stopCount = new AtomicInteger(count)
        final remainder = opts.remainder ? opts.remainder as boolean : false

        DataflowHelper.subscribeImpl( source, phaseHandler(state, count, 0, result, mapper, stopCount, remainder) )
        DataflowHelper.subscribeImpl( target, phaseHandler(state, count, 1, result, mapper, stopCount, remainder) )
        return result
    }

    /**
     * Returns the methods {@code OnNext} and {@code onComplete} which will implement the phase logic
     *
     * @param buffer The shared state buffering the channel received values
     * @param count The overall number of channel
     * @param current The current channel
     * @param target The channel over which the results are sent
     * @param mapper A closure mapping a value to its key
     * @return A map with {@code OnNext} and {@code onComplete} methods entries
     */
    static private final Map phaseHandler( Map<Object,Map<Integer,List>> buffer, int size, int index, DataflowWriteChannel target, Closure mapper, AtomicInteger stopCount, boolean remainder ) {

        [
                onNext: {
                    synchronized (buffer) {
                        def entries = phaseImpl(buffer, size, index, it, mapper, false)
                        if( entries ) {
                            target.bind(entries)
                        }
                    }},

                onComplete: {
                    if( stopCount.decrementAndGet()==0) {
                        if( remainder )
                            phaseRemainder(buffer,size, target)
                        target << Channel.STOP
                    }}

        ]

    }


    /**
     * Implements the phase operator logic. Basically buffers the values received on each channel by their key .
     *
     * When a value with the same key has arrived on each channel, they are removed from the buffer and returned as list
     *
     *
     * @param buffer The shared state buffer
     * @param size The overall number of channels
     * @param current The current channel
     * @param item The value just arrived
     * @param mapper The mapping closure retrieving a key by the item just arrived over the current channel
     * @return The list a values having a common key for each channel or {@code null} when some values are missing
     *
     */
    @PackageScope
    static final List phaseImpl( Map<Object,Map<Integer,List>> buffer, int size, int index, def item, Closure mapper, boolean isCross = false) {

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
        final key = mapper.call(item)

        // given a key we expect to receive on object with the same key on each channel
        def channels = buffer.get(key)
        if( channels==null ) {
            channels = new TreeMap<Integer, List>()
            buffer[key] = channels
        }

        if( !channels.containsKey(index) ) {
            channels[index] = []
        }
        def entries = channels[index]

        // add the received item to the list
        // when it is used in the gather op add always as the first item
        if( isCross && index == 0 ) {
            entries[0] = item
        }
        else  {
            entries << item
        }

        // now check if it has received a element matching for each channel
        if( channels.size() != size )  {
            return null
        }

        def result = []

        Iterator<Map.Entry<Integer,List>> itr = channels.iterator()
        while( itr.hasNext() ) {
            def entry = itr.next()

            def list = entry.getValue()
            result << list[0]

            // do not remove the first element when it is 'cross' op
            if( isCross && entry.getKey() == 0 )
                continue

            list.remove(0)
            if( list.size() == 0 ) {
                itr.remove()
            }
        }

        return result
    }


    static private final void phaseRemainder( Map<Object,Map<Integer,List>> buffers, int count, DataflowWriteChannel target ) {
        Collection<Map<Integer,List>> slots = buffers.values()

        slots.each { Map<Integer,List> entry ->

            while( true ) {

                boolean fill=false
                def result = new ArrayList(count)
                for( int i=0; i<count; i++ ) {
                    List values = entry[i]
                    if( values ) {
                        fill |= true
                        result[i] = values[0]
                        values.remove(0)
                    }
                    else {
                        result[i] = null
                    }
                }

                if( fill )
                    target.bind( result )
                else
                    break
            }

        }
    }

}
