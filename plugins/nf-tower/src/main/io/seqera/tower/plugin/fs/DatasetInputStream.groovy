/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin.fs

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel

import groovy.transform.CompileStatic

/**
 * Minimal {@link SeekableByteChannel} backed by an {@link InputStream}.
 * Supports sequential reads only (no seek/position).
 *
 * @author Seqera Labs
 */
@CompileStatic
class DatasetInputStream implements SeekableByteChannel {
    private final InputStream inputStream
    private long position0 = 0L
    private boolean open = true
    private byte[] buf = new byte[0]

    DatasetInputStream(InputStream inputStream) {
        this.inputStream = inputStream
    }

    @Override
    int read(ByteBuffer dst) throws IOException {
        final len = dst.remaining()
        if (buf.length < len)
            buf = new byte[len]
        final n = inputStream.read(buf, 0, len)
        if (n > 0) {
            dst.put(buf, 0, n)
            position0 += n
        }
        return n
    }

    @Override
    int write(ByteBuffer src) { throw new UnsupportedOperationException() }

    @Override
    long position() { position0 }

    @Override
    SeekableByteChannel position(long newPosition) { throw new UnsupportedOperationException("seek not supported") }

    @Override
    long size() { throw new UnsupportedOperationException("size not available for streaming dataset channel") }

    @Override
    SeekableByteChannel truncate(long size) { throw new UnsupportedOperationException() }

    @Override
    boolean isOpen() { open }

    @Override
    void close() throws IOException {
        open = false
        inputStream.close()
    }
}
