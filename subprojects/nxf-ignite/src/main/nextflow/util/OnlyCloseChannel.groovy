/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
