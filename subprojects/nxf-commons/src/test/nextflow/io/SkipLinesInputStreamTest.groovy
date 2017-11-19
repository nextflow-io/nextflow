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

    def 'should keep the header' () {

        given:
        def header
        def lines
        SkipLinesInputStream filter
        def TEXT = 'FOO\nAA\nBB\rCC\r\nDD\n\rZZ'

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 1)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\n'
        lines == ['AA','BB','CC','DD','','ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 2)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\n'
        lines == ['BB','CC','DD','','ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 3)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\r'
        lines == ['CC','DD','','ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 4)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\rCC\r\n'
        lines == ['DD','','ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 5)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\rCC\r\nDD\n'
        lines == ['','ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 6)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\rCC\r\nDD\n\r'
        lines == ['ZZ']

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 7)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\rCC\r\nDD\n\rZZ'
        lines == []

        when:
        filter = new SkipLinesInputStream(new ByteArrayInputStream(TEXT.bytes), 10)
        header = filter.consumeHeader()
        lines = filter.readLines()
        then:
        header == 'FOO\nAA\nBB\rCC\r\nDD\n\rZZ'
        lines == []

    }

}
