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

/**
 * Hold a buffer for transfer a remote object chunk
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChunkBuffer implements Comparable<ChunkBuffer> {

    private ByteBuffer target;

    private final ChunkBufferFactory owner;

    private volatile int index;

    ChunkBuffer(ChunkBufferFactory owner, int capacity, int index) {
        this.owner = owner;
        this.target = ByteBuffer.allocateDirect(capacity);
        this.index = index;
    }

    private ChunkBuffer(ByteBuffer buffer) {
        this.target = buffer;
        this.owner = null;
    }

    ChunkBuffer withIndex(int index) {
        this.index = index;
        return this;
    }

    int getIndex() {
        return index;
    }

    int getByte() {
        return target.get() & 0xFF;
    }

    void writeByte(int ch) {
        target.put((byte)ch);
    }

    void fill(InputStream stream) throws IOException {
        int ch;
        while ((ch = stream.read())!=-1 && !Thread.currentThread().isInterrupted()) {
            this.writeByte(ch);
        }
    }

    void makeReadable() {
        target.flip();
    }

    void clear() {
        target.clear();
    }

    boolean hasRemaining() {
        return target.hasRemaining();
    }

    public void release() {
        if( owner!=null )
            owner.giveBack(this);

    }

    public static ChunkBuffer wrap(byte[] data) {
        return new ChunkBuffer(ByteBuffer.wrap(data));
    }

    @Override
    public int compareTo(ChunkBuffer other) {
        return Integer.compare(index, other.index);
    }
}
