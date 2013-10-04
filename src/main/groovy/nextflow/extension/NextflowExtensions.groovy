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

package nextflow.extension

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicLong

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
/**
 * Provides extension methods to chunk text and file
 *
 * See more about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class NextflowExtensions {


    def static String rightTrim(String self) {
        self.replaceAll(/\s+$/,"")
    }

    def static String leftTrim( String self ) {
        self.replaceAll(/^\s+/,"")
    }

    /**
     * Splits a {@code CharSequence} in text chunks having the specified number of lines
     *
     * @param sequence A sequence of chars to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    public static void chunkLines( CharSequence sequence, int n = 1, Closure block ) {
        assert sequence != null
        chunkLines( new StringReader(sequence.toString()), [size: n], block )
    }

    /**
     * Splits a {@code CharSequence} in text chunks having the specified number of lines
     *
     * @param sequence A sequence of chars to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    public static void chunkLines( CharSequence sequence, Map options, Closure block ) {
        assert sequence != null
        chunkLines( new StringReader(sequence.toString()), options, block )
    }

    /**
     * Splits a {@code Reader} in text chunks having the specified number of lines
     *
     * @param reader A reader to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    static void chunkLines( Reader reader, int n = 1,  Closure block) {
        chunkLines( reader, [size: n], block)
    }


    /**
     * Splits a {@code Reader} in text chunks having the specified number of lines
     *
     * @param reader A reader to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */

    static void chunkLines( Reader reader, Map options,  Closure block) {
        assert reader != null
        assert options != null
        assert block != null

        log.debug "Chunk options: ${options}"

        int size = options?.size ?: 1
        log.debug "Chunk size: $size"

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(block)

        String line
        StringBuilder buffer = new StringBuilder()
        int c=0
        while( (line = reader0.readLine()) != null ) {
            if ( c ) buffer << '\n'
            buffer << line
            if ( ++c == size ) {
                c = 0
                block.call( buffer.toString() )

                buffer.setLength(0)
            }
        }

        if ( buffer.size() ) {
            block.call( buffer.toString() )
        }

        // send a poison pill to any written channel
        def close = options.autoClose
        if( close == null || close == true ) {
            interceptor.closeChannels()
        }
        else {
            log.debug "Skipping channel autoClose"
        }

    }


    /**
     * Splits a {@code File} in text chunks having the specified number of lines
     *
     * @param file A file to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    static void chunkLines( File file, int n = 1, Closure block ) {
        assert file
        chunkLines( new FileReader(file), [size: n], block )
    }

    /**
     * Splits a {@code File} in text chunks having the specified number of lines
     *
     * @param file A file to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    static void chunkLines( File file, Map options, Closure block ) {
        assert file
        chunkLines( new FileReader(file), options, block )
    }

    static void chunkLines( Path file, int n = 1, Closure block ) {
        assert file
        chunkLines( Files.newBufferedReader(file, Charset.defaultCharset()), n, block )
    }


    static void chunkLines( Path file, Map options, Closure block ) {
        assert file
        chunkLines( Files.newBufferedReader(file, Charset.defaultCharset()), options, block )
    }

    /**
     * Splits a {@code InputStream} in text chunks having the specified number of lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    static void chunkLines( InputStream stream, int n = 1, Closure block ) {
        assert stream
        chunkLines( new InputStreamReader(stream), [size: n], block )
    }

    /**
     * Splits a {@code InputStream} in text chunks having the specified number of lines
     *
     * @param stream A text stream to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */

    static void chunkLines( InputStream stream, Map options, Closure block ) {
        assert stream
        chunkLines( new InputStreamReader(stream), options, block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code CharSequence} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    public static void chunkFasta( CharSequence text, int n = 1, Closure block ) {
        assert text != null
        chunkFasta( new StringReader(text.toString()), [size: n], block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code CharSequence} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */

    public static void chunkFasta( CharSequence text, Map options, Closure block ) {
        assert text != null
        chunkFasta( new StringReader(text.toString()), options, block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code InputStream} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    static void chunkFasta( InputStream text, int n = 1, Closure block ) {
        assert text != null
        chunkFasta( new InputStreamReader(text), [size: n], block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code InputStream} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */

    static void chunkFasta( InputStream text, Map options, Closure block ) {
        assert text != null
        chunkFasta( new InputStreamReader(text), options, block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code File} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    static void chunkFasta( File file, int n = 1, Closure block ) {
        assert file
        chunkFasta( new FileReader(file), n, block )
    }

    static void chunkFasta( Path path, int n = 1, Closure block ) {
        assert path
        chunkFasta( Files.newBufferedReader(path, Charset.defaultCharset()), n, block )
    }


    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code File} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */

    static void chunkFasta( File file, Map options, Closure block ) {
        assert file
        chunkFasta( new FileReader(file), options, block )
    }

    static void chunkFasta( Path file, Map options, Closure block ) {
        assert file
        chunkFasta( Files.newBufferedReader(file, Charset.defaultCharset()), options, block )
    }


    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code Reader} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    static void chunkFasta( Reader text, int n = 1,  Closure block ) {
        chunkFasta( text, [size: n], block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code Reader} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */

    static void chunkFasta( Reader text, Map options, Closure block ) {
        assert text != null
        assert block != null
        assert options != null

        log.debug "Chunk options: $options"

        int size = options.size ?: 1
        log.debug "Chunk size: $size"

        BufferedReader reader0 = text instanceof BufferedReader ? text : new BufferedReader(text)
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

                if ( ++blockCount == size ) {
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
        def close = options.autoClose
        if( close == null || close == true ) {
            interceptor.closeChannels()
        }
        else {
            log.debug "Skipping channel autoClose"
        }
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
                NextflowExtensions.log.error("${e.getMessage()} -- See the file '.nextflow.log' for more error details", e)
                return true
            }
        }

        group.operator(inputs:[queue], outputs:[], listeners: [listenForPoisonPill], code )
    }

    /**
     * Implements a semantically equivalent 'each' iterator over a {@code DataflowStreamReadAdapter} channel

     * @param channel
     * @param closure
     */
    static void eachWithIndex( DataflowStreamReadAdapter channel, Closure closure ) {
        eachWithIndex(channel, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void eachWithIndex( DataflowStreamReadAdapter channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    static void eachWithIndex( DataflowStreamReadAdapter channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }


    static void eachWithIndex( DataflowQueue queue, Closure closure ) {
        eachWithIndex(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void eachWithIndex( DataflowQueue channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    static void eachWithIndex( DataflowQueue channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }


    static void eachWithIndex( WriteChannelWrap channel, Closure code ) {
        eachWithIndex0( channel.target as DataflowReadChannel, Dataflow.retrieveCurrentDFPGroup(), code )
    }


    private static void eachWithIndex0( DataflowReadChannel channel, final PGroup group, Closure code ) {

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(code)

        // -- when a 'PoisonPill' is received, spread it over over any written channel
        def listenForPoisonPill = new DataflowEventAdapter() {
            public Object controlMessageArrived(final DataflowProcessor arg0, final DataflowReadChannel<Object> arg1, final int arg2, final Object message) {
                if ( message instanceof PoisonPill ) {
                    interceptor.getWrittenChannels() *. bind(message)
                }
                return message;
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                NextflowExtensions.log.error("${e.getMessage()} -- See the file '.nextflow.log' for more error details", e)
                return true
            }
        }

        def index = new AtomicLong()
        group.operator(inputs:[channel], outputs:[], listeners: [listenForPoisonPill]) { entry ->
             // invoke the user 'code' passing the current index as an extra parameter
            code.call(entry,index.getAndIncrement())

        }
    }


    /**
     * Keep trace of any reference to a {@code DataflowQueue} or {@code DataflowBroadcast}
     */
    private static class WritableChannelInterceptor {

        /* The target owner object */
        def target

        /* Any reference to a {@code DataflowQueue} or {@code DataflowBroadcast} in the closure block */
        def List<WriteChannelWrap> channels = []

        def private added = []

        def private boolean closed = false

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

        def boolean isClosed() { closed }

        def void closeChannels( ) {
            if( !closed ) {
                channels *.unwrap() .each { channel -> channel << PoisonPill.instance }
                closed = true
            }
        }

        def getProperty(String name) {

            def result = target.getProperty(name)
            if( result instanceof DataflowQueue && !added.contains(result)) {
                added << result
                channels << ( result = new WriteChannelWrap(result) )
            }
            else if ( result instanceof DataflowBroadcast && !added.contains(result) ) {
                added << result
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


}
