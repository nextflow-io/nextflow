/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.util

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ArrayBagTest extends Specification {

    def testOrder( ) {

        given:
        def alpha = new ArrayBag(['abc',123,'x', 9])
        def delta = new ArrayBag([123,9,'abc','x'])

        expect:
        //alpha == delta
        //alpha.hashCode() == delta.hashCode()
        CacheHelper.hasher(alpha).hash() == CacheHelper.hasher(delta).hash()
        CacheHelper.hasher([1,alpha]).hash() == CacheHelper.hasher([1,delta]).hash()

    }

    def testGetAt() {

        given:
        def alpha = new ArrayBag(['abc',123,'x', 9])

        expect:
        alpha[0] == 'abc'
        alpha[1] == 123
        alpha[2] == 'x'
        alpha[3] == 9

    }

    def testSetAt() {

        given:
        def alpha = new ArrayBag(['abc',123,'x', 9])

        when:
        alpha[2] = 'Hello'
        then:
        alpha[0] == 'abc'
        alpha[1] == 123
        alpha[2] == 'Hello'
        alpha[3] == 9

    }

    def 'should stringify the bag' () {
        given:
        def bag = new ArrayBag([1,2,3])
        
        expect:
        bag.toString() == '[1, 2, 3]'
        String.valueOf(bag) == '[1, 2, 3]'
    }

    def 'hashCode should be invariant to order' () {
        given:
        def bag1 = new ArrayBag([1,2,3])
        def bag2 = new ArrayBag([3,1,2])
        def bag3 = new ArrayBag([4,1,2])

        expect:
        bag1.hashCode() == bag2.hashCode()
        bag1.hashCode() != bag3.hashCode()

        /**
         * NOTE!!! equality cannot be checked due to groovy overriding the equals implementation
         * see {@link ArrayBag#equals(java.lang.Object)}
         */
    }

    def 'should access map entry using bag as key' () {
        given:
        def bag1 = new ArrayBag([1,2,3])
        def bag2 = new ArrayBag([3,1,2])
        def bag3 = new ArrayBag([4,1,2])
        and:
        def map = [(bag1):'foo']

        expect:
        map.get(bag1) == 'foo'
        map.get(bag2) == 'foo'
        map.get(bag3) == null

    }

}
