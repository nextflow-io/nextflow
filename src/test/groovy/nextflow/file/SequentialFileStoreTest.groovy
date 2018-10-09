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

import org.junit.Rule
import spock.lang.Specification
import test.TemporaryPath
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SequentialFileStoreTest extends Specification {

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def 'should write and read the same data' () {

        given:
        Path file = tempDir.newFile()
        def HELLO = 'Hello world'.getBytes()

        when:
        def store = new SequentialFileStore(file)
        store.writeBool(true)
        store.writeChar('x' as char)
        store.writeShort( 10 as short )
        store.writeInt( 20 )
        store.writeLong( 30 )
        store.writeBytes( HELLO )
        store.flip()

        then:
        store.readBool()
        store.readChar() == 'x' as char
        store.readShort() == 10 as short
        store.readInt() == 20
        store.readLong() == 30

        def bytes = new byte[HELLO.size()]
        store.readBytes(bytes)
        new String(bytes) == 'Hello world'

        when:
        store.readByte()
        then:
        thrown(EOFException)

        cleanup:
        store.close()

    }

    def 'should write and read the same more data' () {

        given:
        Path file = tempDir.newFile()
        def HELLO = 'Hello world'.getBytes()
        def store = new SequentialFileStore(file, 10) // <- use a very small buffer, this forces to flush the data to the file

        when:
        10.times {

            store.writeBool(true)
            store.writeChar('x' as char)
            store.writeShort( 10 as short )
            store.writeInt( 20 )
            store.writeLong( 30 )
            store.writeBytes( HELLO )

        }
        store.flip()

        then:
        10.times {

            assert store.readBool()
            assert store.readChar() == 'x' as char
            assert store.readShort() == 10 as short
            assert store.readInt() == 20
            assert store.readLong() == 30

            def bytes = new byte[HELLO.size()]
            store.readBytes(bytes)
            assert  new String(bytes) == 'Hello world'

        }

        // check the file has been written and has the same size of written data
        file.size() == store.size()

        cleanup:
        store.close()

    }


    def 'should read and allocate the result byte array' () {

        given:
        Path file = tempDir.newFile()
        def HELLO = 'Hello world'

        when:
        def store = new SequentialFileStore(file)
        store.writeBytes( HELLO.bytes )
        store.writeBytes( HELLO.reverse().bytes )
        store.flip()

        then:
        new String(store.readBytes(HELLO.bytes.size())) == HELLO
        new String(store.readBytes(HELLO.bytes.size())) == HELLO.reverse()

        cleanup:
        store.close()

    }

}
