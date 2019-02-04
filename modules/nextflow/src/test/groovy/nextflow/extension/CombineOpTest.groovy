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
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Session
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CombineOpTest extends Specification {

    def setupSpec() {
        new Session()
    }


    def 'should make a tuple' () {
        given:
        def op = new CombineOp(Mock(DataflowQueue), Mock(DataflowQueue))

        expect:
        op.tuple(pivot, left,right) == result

        where:
        pivot       | left              | right             | result
        [1]         | 'a'               | 'p'               | [1, 'a','p']
        [1]         | ['a','b']         | 'p'               | [1, 'a','b','p']
        [1]         | ['a','b']         | ['p']             | [1, 'a','b','p']
        [1]         | ['a','b']         | []                | [1, 'a','b']
        [1]         | 'a'               | ['p','q']         | [1, 'a', 'p','q']
        [1]         | ['a']             | ['p','q']         | [1, 'a', 'p','q']
        [1]         | []                | ['p','q']         | [1, 'p','q']
        [1]         | ['a','b']         | ['p','q']         | [1, 'a','b','p','q']
        [1]         | ['a',['b','c']]   | 'z'               | [1, 'a', ['b','c'], 'z']
        [2]         | 'a'               | ['p',['w','z']]   | [2, 'a', 'p', ['w','z']]
        [1,2]       | ['a','b']         | ['p','q']         | [1,2, 'a','b','p','q']
    }

    @Timeout(1)
    def 'should combine channels' () {

        given:
        def left = Channel.from('a','b','c')
        def right = Channel.from(1,2,3)
        def op = new CombineOp(left, right)

        when:
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 9
        ['a', 1] in all
        ['a', 2] in all
        ['a', 3] in all
        ['b', 1] in all
        ['b', 2] in all
        ['b', 3] in all
        ['c', 1] in all
        ['c', 2] in all
        ['c', 3] in all
    }

    def 'should combine by a value' () {
        given:
        def left = Channel.from(['A', 10], ['A', 20], ['B', 30], ['B', 40])
        def right = Channel.from(['A', 11], ['A', 22], ['B', 33])

        when:
        def op = new CombineOp(left,right)
        op.setPivot([0])
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val
        println all
        then:
        true
    }

    @Timeout(1)
    def 'should combine a channel with a list' () {

        given:
        def left = Channel.from('a','b')
        def right = [1,2,3,4]
        def op = new CombineOp(left, right)

        when:
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 8
        ['a', 1] in all
        ['a', 2] in all
        ['a', 3] in all
        ['a', 4] in all
        ['b', 1] in all
        ['b', 2] in all
        ['b', 3] in all
        ['b', 4] in all
    }

    @Timeout(1)
    def 'should combine a value with a list' () {

        given:
        def left = Channel.value('x')
        def right = [1,2,3,4]
        def op = new CombineOp(left, right)

        when:
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 4
        ['x', 1] in all
        ['x', 2] in all
        ['x', 3] in all
        ['x', 4] in all
    }

    @Timeout(1)
    def 'should combine two values' () {

        given:
        def left = Channel.value('x')
        def right = Channel.value('z')
        def op = new CombineOp(left, right)

        when:
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 1
        ['x', 'z'] in all
    }

    @Timeout(1)
    def 'should combine with empty value' () {
        given:
        def left = Channel.value().close()
        def right = Channel.value('z')
        def op = new CombineOp(left, right)
        when:
        def result = op.apply()
        then:
        result.val == Channel.STOP
    }

    @Timeout(1)
    def 'should chain combine ops flat default' () {

        given:
        def ch2 = Channel.from('a','b','c')
        def ch3 = Channel.from('x','y')
        def ch1 = Channel.from(1,2)

        when:
        def result = new CombineOp(new CombineOp(ch1, ch2).apply(), ch3).apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 12

        [1,'a','x'] in all
        [1,'a','y'] in all
        [1,'b','x'] in all
        [1,'b','y'] in all
        [1,'c','x'] in all
        [1,'c','y'] in all

        [2,'a','x'] in all
        [2,'a','y'] in all
        [2,'b','x'] in all
        [2,'b','y'] in all
        [2,'c','x'] in all
        [2,'c','y'] in all

    }

    def 'should chain combine ops flat true' () {

        given:
        def ch1 = Channel.from(1,2)
        def ch2 = Channel.from('a','b','c')
        def ch3 = Channel.from('x','y')

        when:
        def result = new CombineOp(new CombineOp(ch1, ch2).apply(), ch3).apply()
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 12

        [1,'a','x'] in all
        [1,'a','y'] in all
        [1,'b','x'] in all
        [1,'b','y'] in all
        [1,'c','x'] in all
        [1,'c','y'] in all

        [2,'a','x'] in all
        [2,'a','y'] in all
        [2,'b','x'] in all
        [2,'b','y'] in all
        [2,'c','x'] in all
        [2,'c','y'] in all
    }

    def 'should combine with tuples' () {

        when:
        def left = Channel.from([1, 'x'], [2,'y'], [3, 'z'])
        def right = ['alpha','beta','gamma']

        def result = new CombineOp(left, right).apply()
        def all = (List) ToListOp.apply(result).val

        then:
        all.size() == 9

        [1, 'x', 'alpha'] in all
        [1, 'x', 'beta'] in all
        [1, 'x', 'gamma'] in all

        [2, 'y', 'alpha'] in all
        [2, 'y', 'beta'] in all
        [2, 'y', 'gamma'] in all

        [3, 'z', 'alpha'] in all
        [3, 'z', 'beta'] in all
        [3, 'z', 'gamma'] in all

    }

    def 'should combine with map' () {

        when:
        def left = Channel.from([id:1, val:'x'], [id:2,val:'y'], [id:3, val:'z'])
        def right = ['alpha','beta','gamma']
        def result = left.combine(right)
        def all = (List) ToListOp.apply(result).val

        then:
        all.size() == 9
        [[id:1, val:'x'], 'alpha'] in all
        [[id:1, val:'x'], 'beta'] in all
        [[id:1, val:'x'], 'gamma'] in all

        [[id:2,val:'y'], 'alpha'] in all
        [[id:2,val:'y'], 'beta'] in all
        [[id:2,val:'y'], 'gamma'] in all

        [[id:3, val:'z'], 'alpha'] in all
        [[id:3, val:'z'], 'beta'] in all
        [[id:3, val:'z'], 'gamma'] in all

    }



    def 'should combine items'() {

        when:
        def left = Channel.from(1,2,3)
        def right = ['a','b']
        def result = left.combine(right).toSortedList().val.iterator()
        then:
        result.next() == [1, 'a']
        result.next() == [1, 'b']
        result.next() == [2, 'a']
        result.next() == [2, 'b']
        result.next() == [3, 'a']
        result.next() == [3, 'b']

        when:
        left = Channel.from(1,2)
        right = Channel.from('a','b','c')
        result = left.combine(right).toSortedList().val.iterator()
        then:
        result.next() == [1, 'a']
        result.next() == [1, 'b']
        result.next() == [1, 'c']
        result.next() == [2, 'a']
        result.next() == [2, 'b']
        result.next() == [2, 'c']

    }

    def 'should chain combine'() {

        when:
        def str1 = Channel.from('a','b','c')
        def str2 = Channel.from('x','y')
        def result = Channel.from(1,2).combine(str1).combine(str2).toSortedList().val.iterator()
        then:
        result.next() == [1,'a','x']
        result.next() == [1,'a','y']
        result.next() == [1,'b','x']
        result.next() == [1,'b','y']
        result.next() == [1,'c','x']
        result.next() == [1,'c','y']
        result.next() == [2,'a','x']
        result.next() == [2,'a','y']
        result.next() == [2,'b','x']
        result.next() == [2,'b','y']
        result.next() == [2,'c','x']
        result.next() == [2,'c','y']

        when:
        str1 = Channel.from('a','b','c')
        str2 = Channel.from('x','y')
        result = Channel.from(1,2).combine(str1).combine(str2,flat:false).toSortedList().val.iterator()
        then:
        result.next() == [1,'a','x']
        result.next() == [1,'a','y']
        result.next() == [1,'b','x']
        result.next() == [1,'b','y']
        result.next() == [1,'c','x']
        result.next() == [1,'c','y']
        result.next() == [2,'a','x']
        result.next() == [2,'a','y']
        result.next() == [2,'b','x']
        result.next() == [2,'b','y']
        result.next() == [2,'c','x']
        result.next() == [2,'c','y']
    }

    def 'should combine by first element' () {

        given:
        def left = Channel.from( ['A',1], ['A',2], ['B',1], ['B',2] )
        def right = Channel.from( ['A',1], ['A',2], ['B',1], ['B',2] )

        when:
        def op = new CombineOp(left, right)
        op.pivot = 0
        def result = op.apply()
        def all = (List) ToListOp.apply(result).val

        then:
        all.size() == 8

        ['A', 1,1] in all
        ['A', 1,2] in all
        ['A', 2,1] in all
        ['A', 2,2] in all

        ['B', 1,1] in all
        ['B', 1,2] in all
        ['B', 2,1] in all
        ['B', 2,2] in all
    }

    def 'should combine by first item' () {
        given:
        def left = [[1, 'a'], [2,'b'], [1,'c']]
        def right = [[2,'p'], [2,'q'], [1,'r']]

        when:
        def result = left.channel().combine(right.channel())
        def all = (List) ToListOp.apply(result).val
        then:
        all.size() == 9
        [1, 'a', 2, 'p'] in all
        [1, 'a', 2, 'q'] in all
        [1, 'a', 1, 'r'] in all

        [2, 'b', 2, 'p'] in all
        [2, 'b', 2, 'q'] in all
        [2, 'b', 1, 'r'] in all

        [1, 'c', 2, 'p'] in all
        [1, 'c', 2, 'q'] in all
        [1, 'c', 1, 'r'] in all

        when:
        result = left.channel().combine(right.channel(), by: 0)
        all = (List) ToListOp.apply(result).val
        then:
        all.size() == 4
        [1, 'a', 'r'] in all
        [1, 'c', 'r'] in all
        [2, 'b', 'p'] in all
        [2, 'b', 'q'] in all

        when:
        result = left.channel().combine(right.channel(), by: [0])
        all = (List) ToListOp.apply(result).val
        then:
        all.size() == 4
        [1, 'a', 'r'] in all
        [1, 'c', 'r'] in all
        [2, 'b', 'p'] in all
        [2, 'b', 'q'] in all
    }

}
