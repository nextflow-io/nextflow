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
 */

package nextflow.io
import java.nio.ByteBuffer

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ByteBufferBackedInputStreamTest extends Specification {

    def testAdaptor () {

        setup:
        def buffer = ByteBuffer.wrap( "Hello\nworld!".getBytes() )

        when:
        def stream = new BufferedInputStream(new ByteBufferBackedInputStream(buffer))
        def lines = stream.readLines()

        then:
        lines[0] == 'Hello'
        lines[1] == 'world!'
        lines.size() == 2

    }

}
