/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

    def 'should join a singleton value' () {

        when:
        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.value(1)

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result == [1]
    }



}
