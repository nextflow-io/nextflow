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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ArrayTupleTest extends Specification {

    def testTuple() {
        when:
        def tuple = new ArrayTuple([1,2,3])
        then:
        tuple[0] == 1
        tuple[1] == 2
        tuple[2] == 3
        tuple.size() == 3

        when:
        tuple.add(4)
        then:
        thrown(UnsupportedOperationException)

        when:
        tuple.remove(0)
        then:
        thrown(UnsupportedOperationException)

        when:
        tuple[0] = 2
        then:
        thrown(UnsupportedOperationException)

        expect:
        tuple.sum() == 6
        tuple.iterator().next() == 1


    }
}
