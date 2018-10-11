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

import spock.lang.Timeout

import nextflow.Channel
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class BufferOpTest extends Specification {

    
    def testBufferClose() {

        when:
        def r1 = Channel.from(1,2,3,1,2,3).buffer({ it == 2 })
        then:
        r1.val == [1,2]
        r1.val == [3,1,2]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from('a','b','c','a','b','z').buffer(~/b/)
        then:
        r2.val == ['a','b']
        r2.val == ['c','a','b']
        r2.val == Channel.STOP

    }

    def testBufferWithCount() {

        when:
        def r1 = Channel.from(1,2,3,1,2,3,1).buffer( size:2 )
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == Channel.STOP

        when:
        r1 = Channel.from(1,2,3,1,2,3,1).buffer( size:2, remainder: true )
        then:
        r1.val == [1,2]
        r1.val == [3,1]
        r1.val == [2,3]
        r1.val == [1]
        r1.val == Channel.STOP


        when:
        def r2 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2,9).buffer( size:3, skip:2 )
        then:
        r2.val == [3,4,5]
        r2.val == [3,4,5]
        r2.val == Channel.STOP

        when:
        r2 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2,9).buffer( size:3, skip:2, remainder: true )
        then:
        r2.val == [3,4,5]
        r2.val == [3,4,5]
        r2.val == [9]
        r2.val == Channel.STOP

    }

    def testBufferInvalidArg() {

        when:
        Channel.create().buffer( xxx: true )

        then:
        IllegalArgumentException e = thrown()

    }


    def testBufferOpenClose() {

        when:
        def r1 = Channel.from(1,2,3,4,5,1,2,3,4,5,1,2).buffer( 2, 4 )
        then:
        r1.val == [2,3,4]
        r1.val == [2,3,4]
        r1.val == Channel.STOP

        when:
        def r2 = Channel.from('a','b','c','a','b','z').buffer(~/a/,~/b/)
        then:
        r2.val == ['a','b']
        r2.val == ['a','b']
        r2.val == Channel.STOP

    }

    def testBufferCloseWithOptions() {

        when:
        def sum = 0
        def r1 = Channel.from(1,2,3,1,2,3).buffer(remainder: true, { sum+=it; sum==7 })
        then:
        r1.val == [1,2,3,1]
        r1.val == [2,3]
        r1.val == Channel.STOP

    }

    def testBufferWithValueChannel() {

        when:
        def result = Channel.value(1).buffer(size: 1)
        then:
        result.val == [1]
        result.val == Channel.STOP

        when:
        result = Channel.value(1).buffer(size: 10)
        then:
        result.val == Channel.STOP

        when:
        result = Channel.value(1).buffer(size: 10,remainder: true)
        result.val == [1]
        then:
        result.val == Channel.STOP
    }


}
