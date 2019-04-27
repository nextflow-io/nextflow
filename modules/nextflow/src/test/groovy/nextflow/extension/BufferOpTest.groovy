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

package nextflow.extension

import spock.lang.Timeout

import nextflow.Channel
import spock.lang.Specification

import nextflow.Session

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class BufferOpTest extends Specification {

    def setup() {
        new Session()
    }
    
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
