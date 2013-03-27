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

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
/**
 * Provides extension methods to chunk text and file
 *
 * See about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChunkMethodsExtension {

    /**
     * Splits a {@code CharSequence} by lines
     *
     * @param sequence The sequence of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    public static DataflowWriteChannel chunkLines( CharSequence sequence, int n = 1, DataflowWriteChannel channel = new DataflowQueue() ) {
        assert sequence
        chunkLines( new StringReader(sequence.toString()), n, channel )
    }

    /**
     * Splits the content of {@code Reader} by lines
     *
     * @param reader The reader of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkLines( Reader reader, int n = 1,  DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert reader

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

        String line
        StringBuilder buffer = new StringBuilder()
        int c=0
        while( (line = reader0.readLine()) != null ) {
            if ( c ) buffer << '\n'
            buffer << line
            if ( ++c == n ) {
                c = 0
                channel << buffer.toString()
                buffer.setLength(0)
            }
        }

        if ( buffer.size() ) {
            channel << buffer.toString()
        }

        // send a poison pill to 'close' the channel
        channel << PoisonPill.instance

        channel
    }


    /**
     * Splits the content of {@code InputStream} by lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkLines( File file, int n = 1, DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert file
        chunkLines( new FileReader(file), n, channel )
    }


    /**
     * Splits the content of {@code InputStream} by lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkLines( InputStream stream, int n = 1, DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert stream
        chunkLines( new InputStreamReader(stream), n, channel)
    }

    /**
     * Split a sequence of char FASTA formatted in chunks
     *
     * @param sequence A {@code CharSequence} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    public static DataflowWriteChannel chunkFasta( CharSequence sequence, int n = 1, DataflowWriteChannel channel = new DataflowQueue() ) {
        assert sequence
        chunkFasta( new StringReader(sequence.toString()), n, channel )
    }

    /**
     * Split a {@code InputStream} containing FASTA formatted data in chunks
     *
     * @param sequence A {@code InputStream} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkFasta( InputStream stream, int n = 1, DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert stream
        chunkFasta( new InputStreamReader(stream), n, channel)
    }

    /**
     * Split a {@code File} containing FASTA formatted data in chunks
     *
     * @param sequence A {@code InputStream} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkFasta( File file, int n = 1, DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert file
        chunkFasta( new FileReader(file), n, channel )
    }

    /**
     * Split a {@code InputStream} containing FASTA formatted data in chunks
     *
     * @param sequence A {@code InputStream} to which apply the operation
     * @param n The number of (FASTA) sequences in each chunk
     * @param channel The channel to which send the result, use {@code null} to create a new {@code DataflowQueue} channel instance
     * @return The channel to which chunks are written
     */
    static DataflowWriteChannel chunkFasta( Reader reader, int n = 1,  DataflowWriteChannel channel = new DataflowQueue()  ) {
        assert reader

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

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
                    channel << buffer.toString()
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
            channel << buffer.toString()
        }

        // send a poison pill to close the channel
        channel << PoisonPill.instance

        channel
    }





}
