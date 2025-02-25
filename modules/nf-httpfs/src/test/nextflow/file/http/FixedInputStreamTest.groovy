/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.file.http

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FixedInputStreamTest extends Specification {

    def 'should read byte by byte' () {
        given:
        def bytes = "Hello world". bytes
        def stream = new FixedInputStream(new ByteArrayInputStream(bytes), bytes.length)

        when:
        def ch
        def result = new StringBuilder()
        while( (ch=stream.read())!=-1 )
            result.append(ch as char)
        and:
        stream.close()
        then:
        noExceptionThrown()
        result.toString() == 'Hello world'
    }

    def 'should read byte buffer' () {
        given:
        def bytes = "Hello world". bytes
        def stream = new FixedInputStream(new ByteArrayInputStream(bytes), bytes.length)

        when:
        def buffer = new byte[5]
        def result = new StringBuilder()
        def c
        while( (c=stream.read(buffer))!=-1 ) {
            for( int i=0; i<c; i++ )
                result.append(buffer[i] as char)
        }
        and:
        stream.close()
        then:
        noExceptionThrown()
        result.toString() == 'Hello world'
    }

    def 'should read all bytes' () {
        given:
        def bytes = "Hello world". bytes
        def stream = new FixedInputStream(new ByteArrayInputStream(bytes), bytes.length)

        when:
        def result = stream.readAllBytes()
        and:
        stream.close()
        then:
        noExceptionThrown()
        and:
        new String(result) == "Hello world"
    }
}
