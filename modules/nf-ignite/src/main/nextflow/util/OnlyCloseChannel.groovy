/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel

import groovy.transform.CompileStatic
/**
 * Fake channel only supporting the close operation.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class OnlyCloseChannel implements SeekableByteChannel {

    Closeable channel

    OnlyCloseChannel( Closeable channel) {
        this.channel = channel
    }

    @Override
    int read(ByteBuffer dst) throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    int write(ByteBuffer src) throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    long position() throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    SeekableByteChannel position(long newPosition) throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    long size() throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    SeekableByteChannel truncate(long size) throws IOException {
        throw new UnsupportedOperationException()
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
