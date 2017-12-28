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
