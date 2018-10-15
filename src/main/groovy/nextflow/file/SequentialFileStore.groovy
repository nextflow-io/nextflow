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

package nextflow.file

import java.nio.file.Path

import groovy.transform.CompileStatic
/**
 * A buffered storage. When data written is greater than the local buffer the content is save to a file.
 *
 * <p>
 * Usage idiom:
 *
 * <code>
 * def store = new SequentialFileStore(temp_file)
 * store.writeInt(x)
 * store.writeLong(y)
 * :
 * store.flip()
 * int x = store.readInt()
 * long y = store.readLong()
 * :
 * store.close()
 *
 * </code>
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SequentialFileStore implements Closeable {

    final private Path path

    private RandomAccessFile store

    SequentialFileStore(Path path) {
        this.path = path
        this.store = new RandomAccessFile(path.toFile(), 'rw')
    }

    SequentialFileStore writeBool( boolean value ) {
        store.writeBoolean(value)
        return this
    }

    SequentialFileStore writeByte( byte value ) {
        store.writeByte(value)
        return this
    }

    SequentialFileStore writeChar( char value ) {
        store.writeChar((int)value)
        return this
    }

    SequentialFileStore writeShort( short value ) {
        store.writeShort(value as short)
        return this
    }

    SequentialFileStore writeInt( int v ) {
        store.writeInt(v)
        return this
    }

    SequentialFileStore writeLong( long v ) {
        store.writeLong(v)
        return this
    }


    SequentialFileStore writeBytes( byte[] bytes ) {
        store.write(bytes)
        return this
    }

    void readBytes( byte[] bytes ) {
        store.read(bytes)
    }

    byte[] readBytes( int len ) {
        def bytes = new byte[len]
        store.read(bytes)
        return bytes
    }

    boolean readBool() {
        store.readBoolean()
    }

    byte readByte() {
        store.readByte()
    }

    char readChar() {
        store.readChar()
    }

    short readShort() {
        store.readShort()
    }

    int readInt() {
        store.readInt()
    }

    long readLong() {
        store.readLong()
    }

    long size() {
        store.length()
    }

    /**
     * Turn the buffer in read mode
     */
    void flip() {
        // flush the content to the file
        store.getFD().sync()
        store.seek(0)
    }

    @Override
    void close() throws IOException {
        store.close()
    }
}
