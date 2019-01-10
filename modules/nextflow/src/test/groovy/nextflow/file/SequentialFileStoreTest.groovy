/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
        def store = new SequentialFileStore(file) 

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
