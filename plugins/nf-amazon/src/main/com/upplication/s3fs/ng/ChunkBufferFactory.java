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

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
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
        this.pool = new LinkedBlockingQueue<>();
        this.count = new AtomicInteger();
    }

    public ChunkBuffer create(int index) throws InterruptedException {
        int it=0;
        while( true ) {
            if( (it+1) % 10 == 0 )
                log.debug("Waiting for a download chunk buffer to become available - Consider to decrease the download workers or increase download buffer capacity");
            ChunkBuffer result = it++==0 ? pool.poll() : pool.poll(1, TimeUnit.SECONDS);
            if( result != null ) {
                result.clear();
                return result.withIndex(index);
            }
            if( count.getAndIncrement() < capacity ) {
                return new ChunkBuffer(this,chunkSize,index);
            }
        }
    }

    void giveBack(ChunkBuffer buffer) {
        pool.offer(buffer);
    }

}
