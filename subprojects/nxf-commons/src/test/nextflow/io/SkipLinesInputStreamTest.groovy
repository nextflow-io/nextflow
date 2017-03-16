package nextflow.io

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SkipLinesInputStreamTest extends Specification {

    @Unroll
    def 'should skip first line: #skip' () {

        given:
        def TEXT = 'HEADER\nA\rB\r\nC\n'

        when:
        def filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), SKIP)
        then:
        filter.readLines() == EXPECT

        where:
        SKIP    | EXPECT
        1       | ['A','B','C']
        2       | ['B','C']
        3       | ['C']
        4       | []
        10      | []

    }

}
