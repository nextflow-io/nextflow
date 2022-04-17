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
import java.nio.channels.SeekableByteChannel
import java.nio.channels.WritableByteChannel

import groovy.transform.CompileStatic

/**
 * Implements a {@link SeekableByteChannel} for a given {@link java.nio.channels.WritableByteChannel}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzWriteableByteChannel implements SeekableByteChannel {

    private long _pos
    private WritableByteChannel writer

    AzWriteableByteChannel(WritableByteChannel channel) {
        this.writer = channel
    }

    @Override
    int read(ByteBuffer dst) throws IOException {
        throw new UnsupportedOperationException("Operation 'read(ByteBuffer)' is not supported by AzWriteableByteChannel")
    }

    @Override
    int write(ByteBuffer src) throws IOException {
        def len = writer.write(src)
        _pos += len
        return len
    }

    @Override
    long position() throws IOException {
        return _pos
    }

    @Override
    SeekableByteChannel position(long newPosition) throws IOException {
        throw new UnsupportedOperationException("Operation 'position(long)' is not supported by AzWriteableByteChannel")
    }

    @Override
    long size() throws IOException {
        return _pos
    }

    @Override
    SeekableByteChannel truncate(long size) throws IOException {
        throw new UnsupportedOperationException("Operation 'truncate(long)' is not supported by AzWriteableByteChannel")
    }

    @Override
    boolean isOpen() {
        return true
    }

    @Override
    void close() throws IOException {
        writer.close()
    }

}
