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

package nextflow.lineage.fs

import groovy.transform.CompileStatic

import java.nio.ByteBuffer
import java.nio.channels.ClosedChannelException
import java.nio.channels.NonWritableChannelException
import java.nio.channels.SeekableByteChannel

/**
 * SeekableByteChannel for metadata results description as a file.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class LinMetadataSeekableByteChannel implements SeekableByteChannel {
    private final ByteBuffer buffer
    private boolean open

    LinMetadataSeekableByteChannel(byte[] bytes){
        this.open = true
        this.buffer = ByteBuffer.wrap(bytes)
    }

    @Override
    int read(ByteBuffer dst) {
        if (!open) throw new ClosedChannelException()
        if (!buffer.hasRemaining()) return -1
        int remaining = Math.min(dst.remaining(), buffer.remaining())
        byte[] temp = new byte[remaining]
        buffer.get(temp)
        dst.put(temp)
        return remaining
    }

    @Override
    int write(ByteBuffer src) { throw new NonWritableChannelException() }

    @Override
    long position() { return buffer.position() }

    @Override
    SeekableByteChannel position(long newPosition) {
        if (newPosition < 0 || newPosition > buffer.limit()) throw new IllegalArgumentException()
        buffer.position((int) newPosition)
        return this
    }

    @Override
    long size() { return buffer.limit() }

    @Override
    SeekableByteChannel truncate(long size) { throw new NonWritableChannelException() }

    @Override
    boolean isOpen() { return open }

    @Override
    void close() { open = false }
}
