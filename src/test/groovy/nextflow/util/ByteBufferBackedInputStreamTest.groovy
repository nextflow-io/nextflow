package nextflow.util

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
