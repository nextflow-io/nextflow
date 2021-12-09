/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package com.upplication.s3fs.ng;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Hold a buffer for transfer a remote object chunk
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChunkBuffer {

    static final private Logger log = LoggerFactory.getLogger(ChunkBuffer.class);

    private ByteBuffer target;
    private TransferRateMeter meter = new TransferRateMeter();

    private ChunkBuffer(ByteBuffer buffer) {
        this.target = buffer;
    }

    public static ChunkBuffer create(int capacity) {
        return new ChunkBuffer(ByteBuffer.allocateDirect(capacity));
    }

    void clear() {
        target.clear();
    }

    int getByte() {
        return target.get() & 0xFF;
    }

    void writeByte(int ch) {
        target.put((byte)ch);
    }

    void fill(InputStream stream) throws IOException {
        int ch;
        while ((ch = stream.read()) != -1) {
            this.writeByte(ch);
            meter.inc(1);
        }
    }

    void flip() {
        target.flip();
    }

    boolean hasRemaining() {
        return target.hasRemaining();
    }

    public void release() {
        // return this buffer to teh buffer pool
        target = null;
    }

    public static ChunkBuffer wrap(byte[] data) {
        return new ChunkBuffer(ByteBuffer.wrap(data));
    }

}
