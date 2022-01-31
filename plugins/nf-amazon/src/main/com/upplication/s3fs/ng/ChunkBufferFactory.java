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

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Model a buffer for download chunk
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChunkBufferFactory {

    final Logger log = LoggerFactory.getLogger(ChunkBufferFactory.class);

    final private BlockingQueue<ChunkBuffer> pool;

    final private AtomicInteger count;

    private final int chunkSize;

    private final int capacity;

    public ChunkBufferFactory(int chunkSize, int capacity) {
        this.chunkSize = chunkSize;
        this.capacity = capacity;
        this.pool = new ArrayBlockingQueue<>(capacity);
        this.count = new AtomicInteger();
    }


    public ChunkBuffer create() throws InterruptedException {
        ChunkBuffer result = pool.poll(100, TimeUnit.MILLISECONDS);
        if( result != null ) {
            result.clear();
            return result;
        }

        // add logistic delay to slow down the allocation of new buffer
        // when the request approach or exceed the max capacity
        final int indx = count.getAndIncrement();
        if( log.isTraceEnabled() )
            log.trace("Creating a new buffer index={}; capacity={}", indx, capacity);
        return new ChunkBuffer(this, chunkSize, indx);
    }

    void giveBack(ChunkBuffer buffer) {
        if( pool.offer(buffer) ) {
            if( log.isTraceEnabled() )
                log.trace("Returning buffer {} to pool size={}", buffer.getIndex(), pool.size());
        }
        else {
            int cc = count.decrementAndGet();
            if( log.isTraceEnabled() )
                log.trace("Returning buffer index={} for GC; pool size={}; count={}", buffer.getIndex(), pool.size(), cc);
        }
    }

    int getPoolSize() { return pool.size(); }
}
