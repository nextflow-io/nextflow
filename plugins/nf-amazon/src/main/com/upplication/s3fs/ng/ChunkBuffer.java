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
@Deprecated
public class ChunkBuffer implements Comparable<ChunkBuffer> {

    private static final int BUFFER_SIZE = 8192;

    private final ByteBuffer target;

    private final ChunkBufferFactory owner;

    private final int index;

    ChunkBuffer(ChunkBufferFactory owner, int capacity, int index) {
        this.owner = owner;
        this.target = ByteBuffer.allocateDirect(capacity);
        this.index = index;
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
        int n;
        byte[] b = new byte[BUFFER_SIZE];
        while ((n = stream.read(b)) != -1 ) {
            target.put(b, 0, n);
        }
    }

    void makeReadable() {
        // cast to prevent Java 8 / Java 11 cross compile-runtime error
        // https://www.morling.dev/blog/bytebuffer-and-the-dreaded-nosuchmethoderror/
        ((java.nio.Buffer)target).flip();
    }

    void clear() {
        // cast to prevent Java 8 / Java 11 cross compile-runtime error
        // https://www.morling.dev/blog/bytebuffer-and-the-dreaded-nosuchmethoderror/
        ((java.nio.Buffer)target).clear();
    }

    int getBytes( byte[] buff, int off, int len ) {
        int c=0;
        int i=off;
        while( c<len && target.hasRemaining() ) {
            c++;
            buff[i++] = target.get();
        }
        return c;
    }

    boolean hasRemaining() {
        return target.hasRemaining();
    }

    public void release() {
        owner.giveBack(this);
    }

    @Override
    public int compareTo(ChunkBuffer other) {
        return Integer.compare(index, other.index);
    }
}
