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

package nextflow.extension

import nextflow.Channel
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class JoinOpTest extends Specification {

    def 'should join entries' () {
        given:
        def ch1 = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
        def ch2 = Channel.from(['Z', 6], ['Y', 5], ['X', 4])

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
    }



    def 'should join entries by index' () {
        given:
        def ch1 = Channel.from([1, 'X'], [2, 'Y'], [3, 'Z'], [7, 'P'])
        def ch2 = Channel.from([6, 'Z'], [5, 'Y'], [4, 'X'])

        when:
        def op = new JoinOp(ch1, ch2, [by:1])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
    }

    def 'should join entries with composite index' () {
        given:
        def ch1 = Channel.from([1, 'a','b', ['foo']], [2, 'p','q', ['bar']], [3, 'x','y', ['baz']], [7, 'P'])
        def ch2 = Channel.from([5, 'p','q', [333]], [4, 'a','b', [444]], [6, 'x','y', [555]])

        when:
        def op = new JoinOp(ch1, ch2, [by:[1,2]])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['a','b', 1, ['foo'], 4, [444]] )
        result.contains( ['p','q', 2, ['bar'], 5, [333]] )
        result.contains( ['x','y', 3, ['baz'], 6, [555]] )

    }


    def 'should join entries with remainder' () {
        given:
        def ch1 = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
        def ch2 = Channel.from(['Z', 6], ['Y', 5], ['X', 4], ['Q', ['foo','bar', [77,88,99]]])

        when:
        def op = new JoinOp(ch1, ch2, [remainder: true])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 5
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
        result.contains( ['P', 7, null] )
        result.contains( ['Q', null, ['foo','bar', [77,88,99]]])
    }

    def 'should join single item channels' () {

        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result == [1,2,3]
    }

    def 'should join single item channels with remainder' () {

        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        when:
        def op = new JoinOp(ch1, ch2, [remainder: true])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 8
        result == [1, 2, 3, 0, 0, 7, 8, 9]
    }


}
