/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.nio

import java.nio.ByteBuffer
import java.nio.channels.ReadableByteChannel
import java.nio.channels.SeekableByteChannel

import groovy.transform.CompileStatic

/**
 * Implements a {@link SeekableByteChannel} for a given {@link ReadableByteChannel}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzReadableByteChannel implements SeekableByteChannel {

    private long _position
    private ReadableByteChannel channel
    private long size

    AzReadableByteChannel(ReadableByteChannel channel, long size) {
        this.channel = channel
        this.size = size
    }

    @Override
    int read(ByteBuffer dst) throws IOException {
        final len = channel.read(dst)
        _position += len
        return len
    }

    @Override
    int write(ByteBuffer src) throws IOException {
        throw new UnsupportedOperationException("Operation 'write(ByteBuffer)' is not supported by AzReadableByteChannel")
    }

    @Override
    long position() throws IOException {
        return _position
    }

    @Override
    SeekableByteChannel position(long newPosition) throws IOException {
        throw new UnsupportedOperationException("Operation 'position(long)' is not supported by AzReadableByteChannel")
    }

    @Override
    long size() throws IOException {
        return size
    }

    @Override
    SeekableByteChannel truncate(long dummy) throws IOException {
        throw new UnsupportedOperationException("Operation 'truncate(long)' is not supported by AzReadableByteChannel")
    }

    @Override
    boolean isOpen() {
        return true
    }

    @Override
    void close() throws IOException {
        channel.close()
    }

}
