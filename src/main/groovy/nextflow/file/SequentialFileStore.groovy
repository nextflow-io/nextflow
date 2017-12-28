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
 * A buffered storage. When data written is greater than the local buffer teh content is save to a file.
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

    private static final int KB = 1024

    private static final int DEFAULT_BUFFER_SIZE = 80 * KB

    final private Path path

    private int limit

    private RandomAccessFile file

    private DataOutputStream outputStream

    private DataInputStream inputStream

    private byte[] buffer

    SequentialFileStore(Path path, int size = 0) {
        this.path = path
        this.limit = size ?: DEFAULT_BUFFER_SIZE
        this.buffer = new byte[ size ?: this.limit ]
        this.file = new RandomAccessFile(path.toFile(), 'rw')
        this.outputStream = new DataOutputStream(new FastBufferedOutputStream(new FileOutputStream(file.getFD()), buffer))
    }

    SequentialFileStore writeBool( boolean value ) {
        outputStream.writeBoolean(value)
        return this
    }

    SequentialFileStore writeByte( byte value ) {
        outputStream.writeByte(value)
        return this
    }

    SequentialFileStore writeChar( char value ) {
        outputStream.writeChar(value as int)
        return this
    }

    SequentialFileStore writeShort( short value ) {
        outputStream.writeShort(value as short)
        return this
    }

    SequentialFileStore writeInt( int v ) {
        outputStream.writeInt(v)
        return this
    }

    SequentialFileStore writeLong( long v ) {
        outputStream.writeLong(v)
        return this
    }


    SequentialFileStore writeBytes( byte[] bytes ) {
        outputStream.write(bytes)
        return this
    }

    void readBytes( byte[] bytes ) {
        inputStream.read(bytes)
    }

    byte[] readBytes( int len ) {
        def bytes = new byte[len]
        inputStream.read(bytes)
        return bytes
    }

    boolean readBool() {
        inputStream.readBoolean()
    }

    byte readByte() {
        inputStream.readByte()
    }

    char readChar() {
        inputStream.readChar()
    }

    short readShort() {
        inputStream.readShort()
    }

    int readInt() {
        inputStream.readInt()
    }

    long readLong() {
        inputStream.readLong()
    }

    long size() {
        outputStream.size()
    }

    /**
     * Turn the buffer in read mode
     */
    void flip() {

        if( outputStream.size() > limit ) {
            // flush the content to the file
            outputStream.flush()
            file.seek(0)
            inputStream = new DataInputStream(new FastBufferedInputStream(new FileInputStream(file.getFD()), buffer))
        }
        else {
            // do not flush the buffer and use directly it's content
            inputStream = new DataInputStream( new FastBufferedInputStream(new FastByteArrayInputStream(buffer, 0, outputStream.size())))
        }

    }

    @Override
    void close() throws IOException {
        file.close()
    }
}
