/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package com.upplication.s3fs.ng;

import java.io.IOException;
import java.io.InputStream;
import java.io.InterruptedIOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * Implements an InputStream made of a sequence of chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChunkedInputStream extends InputStream  {

    final private long length;
    private long count;
    private ChunkBuffer buffer;
    private final BlockingQueue<ChunkBuffer> chunks = new LinkedBlockingQueue<>();
    private volatile IOException error;
    private int nextIndex;

    public ChunkedInputStream(long length) {
        this.length = length;
    }

    public void add(byte[] buffer) {
        if( buffer.length>0 ) {
            // skip empty chunks
            chunks.add(ChunkBuffer.wrap(buffer));
        }
    }

    public void add(ChunkBuffer buffer) throws InterruptedException {
        buffer.makeReadable();
        chunks.put(buffer);
    }

    /**
     * Fetch the next bytes from the queue of chunks
     *
     * @return -1 when the end of the stream is reached
     * @throws IOException
     */
    @Override
    public int read() throws IOException {
        if( error != null )
            throw error;
        if( count == length )
            return -1;
        if( buffer == null ) {
            buffer = takeBuffer();
        }
        else if( !buffer.hasRemaining() ) {
            buffer.release();
            buffer = takeBuffer();
        }
        count++;
        return buffer.getByte();
    }

    private ChunkBuffer takeBuffer() throws IOException {
        try {
            while( !Thread.currentThread().isInterrupted() ) {
                if( error != null )
                    throw error;

                ChunkBuffer buffer = chunks.poll(1, TimeUnit.SECONDS);
                if( buffer == null )
                    continue;
                if( buffer.getIndex() != nextIndex ) {
                    chunks.add(buffer);
                    continue;
                }
                nextIndex++;
                return buffer;
            }
            throw new InterruptedIOException("Chunked stream was interrupted");
        }
        catch (InterruptedException e) {
            throw new InterruptedIOException("Chunked stream got interrupted exception");
        }
    }

    public void throwError(IOException t) {
        this.error = t;
    }

    @Override
    public void close() {
        chunks.clear();
        if( buffer!=null )
            buffer.release();
    }
}
