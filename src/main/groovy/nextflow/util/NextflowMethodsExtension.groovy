/*
 * Copyright (c) 2012, the authors.
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

package nextflow.util
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamReadAdapter
import groovyx.gpars.group.DefaultPGroup
import groovyx.gpars.group.PGroup
import groovyx.gpars.scheduler.Pool
import org.apache.commons.io.FileUtils

/**
 * Provides extension methods to chunk text and file
 *
 * See about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class NextflowMethodsExtension {

    /**
     * Splits a {@code CharSequence} by lines
     *
     * @param sequence The sequence of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    public static void chunkLines( CharSequence sequence, int n = 1, Closure block ) {
        assert sequence
        chunkLines( new StringReader(sequence.toString()), n, block )
    }

    /**
     * Splits the content of {@code Reader} by lines
     *
     * @param reader The reader of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static void chunkLines( Reader reader, int n = 1,  Closure block) {
        assert reader

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(block)

        String line
        StringBuilder buffer = new StringBuilder()
        int c=0
        while( (line = reader0.readLine()) != null ) {
            if ( c ) buffer << '\n'
            buffer << line
            if ( ++c == n ) {
                c = 0
                block.call( buffer.toString() )

                buffer.setLength(0)
            }
        }

        if ( buffer.size() ) {
            block.call( buffer.toString() )
        }

        // send a poison pill to any written channel
        interceptor.getWrittenChannels() *.unwrap() .each { channel -> channel << PoisonPill.instance }

    }


    /**
     * Splits the content of {@code InputStream} by lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static void chunkLines( File file, int n = 1, Closure block ) {
        assert file
        chunkLines( new FileReader(file), n, block )
    }


    /**
     * Splits the content of {@code InputStream} by lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static void chunkLines( InputStream stream, int n = 1, Closure block ) {
        assert stream
        chunkLines( new InputStreamReader(stream), n, block )
    }



    public static void chunkFasta( CharSequence sequence, int n = 1, Closure block ) {
        assert sequence
        chunkFasta( new StringReader(sequence.toString()), n, block )
    }

    /**
     * Split a {@code InputStream} containing FASTA formatted data in chunks
     *
     * @param sequence A {@code InputStream} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static void chunkFasta( InputStream stream, int n = 1, Closure block ) {
        assert stream
        chunkFasta( new InputStreamReader(stream), n, block)
    }

    /**
     * Split a {@code File} containing FASTA formatted data in chunks
     *
     * @param sequence A {@code InputStream} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static void chunkFasta( File file, int n = 1, Closure block ) {
        assert file
        chunkFasta( new FileReader(file), n, block )
    }

    static void chunkFasta( Reader reader, int n = 1,  Closure block ) {
        assert reader
        assert block

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)
        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(block)

        String line
        StringBuilder buffer = new StringBuilder()
        int blockCount=0
        boolean openBlock = false
        while( (line = reader0.readLine()) != null ) {

            if ( line == '' ) {
                buffer << '\n'
            }
            else if ( !openBlock && line.charAt(0)=='>' ) {
                openBlock = true
                buffer << line << '\n'
            }
            else if ( openBlock && line.charAt(0)=='>') {
                // another block is started

                if ( ++blockCount == n ) {
                    // invoke the closure, passing the read block as parameter
                    block.call(buffer.toString())

                    buffer.setLength(0)
                    blockCount=0
                }

                buffer << line << '\n'

            }
            else {
                buffer << line << '\n'
            }

        }

        if ( buffer.size() ) {
            block.call(buffer.toString())
        }

        // send a poison pill to any written channel
        interceptor.getWrittenChannels() *.unwrap() .each { channel -> channel << PoisonPill.instance }
    }




    // --==== Dataflow each operator extension ===--


    /**
     * Implements 'each' iterator for {@code DataflowQueue}.
     * <p>
     *     The key feature of this operator is to spread transparently the 'poison pill'
     *     from the source {@code queue} into the referenced queue in the code block
     * <p>
     *     For example:
     *     <code>
     *     source = new DataflowQueue()
     *     target1 = new DataflowQueue()
     *     target2 = new DataflowQueue()
     *
     *     source.each {
     *          target1 << it
     *          target2 << it * it
     *     }
     *     </code>
     *
     *
     *
     * @param queue
     * @param group
     * @param code
     */
    static void each( DataflowStreamReadAdapter queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void each( DataflowStreamReadAdapter queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    static void each( DataflowStreamReadAdapter queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }


    static void each( DataflowQueue queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void each( DataflowQueue queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    static void each( DataflowQueue queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }


    static void each ( WriteChannelWrap channel, Closure code ) {
        each0( channel.target as DataflowReadChannel, Dataflow.retrieveCurrentDFPGroup(), code )
    }

    /**
     * Implements the 'each' operator.
     * <p>
     * The goal is intercept any writes to 'dataflow' queues so that to
     * being able to spread the 'poison pill' from the source 'queue' to the target
     *
     * @param queue A a {@code DataflowQueue} or {@code DataflowBroadcast} object instance
     * @param group
     * @param code The closure code block
     */
    private static void each0( DataflowReadChannel queue, final PGroup group, Closure code ) {

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(code)

        // -- when a 'PoisonPill' is received, spread it over over any written channel
        def listenForPoisonPill = new DataflowEventAdapter() {
            public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if ( message instanceof PoisonPill ) {
                    interceptor.getWrittenChannels() *. bind( message )
                }
                return message;
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                NextflowMethodsExtension.log.error("${e.getMessage()} -- See the file '.nextflow.log' for more error details", e)
                return true
            }
        }

        group.operator(inputs:[queue], outputs:[], listeners: [listenForPoisonPill], code )
    }

    /**
     * Keep trace of any reference to a {@code DataflowQueue} or {@code DataflowBroadcast}
     */
    private static class WritableChannelInterceptor {

        /* The target owner object */
        def target

        /* Any reference to a {@code DataflowQueue} or {@code DataflowBroadcast} in the closure block */
        def List<WriteChannelWrap> channels = []

        WritableChannelInterceptor( Closure code ) {
            assert code
            // replace the closure 'owner' by 'this' instance
            target = code.owner
            code.@owner = this
        }

        /** All the channel for which a 'bind' operation has been invoked */
        List<WriteChannelWrap> getWrittenChannels() {
            channels.findAll{ WriteChannelWrap it -> it.receivedData }
        }

        def getProperty(String name) {

            def result = target.getProperty(name)
            if( result instanceof DataflowQueue ) {
                channels << ( result = new WriteChannelWrap(result) )
            }
            else if ( result instanceof DataflowBroadcast ) {
                channels << ( result = new WriteChannelWrap(result))
            }

            return result
        }

    }

    /**
     * Wrap a {@code WriteChannelWrap} a keep track of any bind operation
     */
    @TupleConstructor
    private static class WriteChannelWrap {

        @Delegate
        DataflowWriteChannel target

        /** Whenever a 'bind' operation has been invoked on the target channel */
        boolean receivedData

        /** The reference to the target channel */
        DataflowWriteChannel unwrap() { target }

        DataflowWriteChannel leftShift(final Object value) {
            receivedData = true
            target.leftShift(value)
        }

        void bind(final Object value) {
            receivedData = true
            target.bind(value)
        }

        DataflowWriteChannel leftShift(final DataflowReadChannel ref) {
            receivedData = true
            target.leftShift(ref)
        }

    }


    // ---=== File helpers ==---


    def static boolean isEmpty( File file ) {
        FileHelper.isEmpty(file)
    }

    def static isNotEmpty( File file ) {
        return FileHelper.isNotEmpty(file)
    }


    def static copyTo( File source, File target ) {
        if( source.isDirectory() ) {
            FileUtils.copyDirectory(source, target)
        }
        else {
            FileUtils.copyFile(source, target)
        }
    }



}
