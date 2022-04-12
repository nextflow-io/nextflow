/*
 * Copyright 2020-2022, Seqera Labs
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
import java.util.Iterator;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Implements an input stream emitting a collection of futures {@link ChunkBuffer}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class FutureInputStream extends InputStream  {

    private final Iterator<Future<ChunkBuffer>> futures;
    private ChunkBuffer buffer;

    FutureInputStream(Iterator<Future<ChunkBuffer>> futures) {
        this.futures = futures;
    }

    @Override
    public int read() throws IOException {

        if( (buffer == null || !buffer.hasRemaining()) ) {
            freeBuffer();
            if( futures.hasNext() )  {
                buffer = nextBuffer();
            }
            else {
                return -1;
            }
        }

        return buffer.getByte();
    }

    @Override
    public int read(byte[] b, int off, int len) throws IOException {

        if( (buffer == null || !buffer.hasRemaining()) ) {
            freeBuffer();
            if( futures.hasNext() )  {
                buffer = nextBuffer();
            }
            else {
                return -1;
            }
        }

        return buffer.getBytes(b, off, len);
    }

    private ChunkBuffer nextBuffer() throws IOException {
        try {
            return futures.next().get();
        }
        catch (ExecutionException e) {
            throw new IOException("Failed to acquire stream chunk", e);
        }
        catch (InterruptedException e) {
            throw new InterruptedIOException();
        }
    }

    private void freeBuffer() {
        if( buffer!=null ) {
            buffer.release();
            buffer=null;
        }
    }

    @Override
    public void close() {
        freeBuffer();
    }
}
