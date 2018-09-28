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
