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
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements an InputStream made of a sequence of chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChunkedInputStream extends InputStream  {

    static final private Logger log = LoggerFactory.getLogger(ChunkedInputStream.class);

    final private long length;
    private long readsCount;
    private ChunkBuffer buffer;
    private final BlockingQueue<ChunkBuffer> chunks = new PriorityBlockingQueue<>();
    private volatile IOException error;
    private int nextIndex;
    private Lock monitor;
    private Condition newChunk;

    public ChunkedInputStream(long length) {
        this.length = length;
        this.monitor = new ReentrantLock();
        this.newChunk = monitor.newCondition();
    }

    public void add(ChunkBuffer chunk) throws InterruptedException {
        monitor.lock();
        try {
            chunks.put(chunk);
            newChunk.signal();
        }
        finally {
            monitor.unlock();
        }
    }

    private void awaitNewChunk() throws InterruptedException {
        monitor.lock();
        try { newChunk.await(10, TimeUnit.MILLISECONDS); }
        finally { monitor.unlock(); }
    }

    /**
     * Fetch the next bytes from the queue of chunks
     *
     * @return -1 when the end of the stream is reached
     * @throws IOException
     */
    @Override
    public int read() throws IOException {
        if( error != null ) {
            throw error;
        }
        if( readsCount == length ) {
            return -1;
        }
        if( buffer == null ) {
            buffer = takeBuffer();
        }
        else if( !buffer.hasRemaining() ) {
            freeBuffer();
            buffer = takeBuffer();
        }
        readsCount++;
        return buffer.getByte();
    }

    private void freeBuffer() {
        if( buffer!=null ) {
            buffer.release();
            buffer=null;
        }
    }

    private ChunkBuffer takeBuffer() throws IOException {
        try {
            while( !Thread.currentThread().isInterrupted() ) {
                if( error != null )
                    throw error;

                ChunkBuffer buffer = chunks.poll(1, TimeUnit.SECONDS);
                if( buffer == null ) {
                    //log.trace("Chunks buffer queue empty - readsCount="+readsCount);
                    continue;
                }
                if( buffer.getIndex() != nextIndex ) {
                    log.trace ("Got chunk with index {}, expected {}", buffer.getIndex(), nextIndex);
                    chunks.add(buffer);
                    awaitNewChunk();
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
        freeBuffer();
    }
}
