/*
 * Copyright 2013-2026, Seqera Labs
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
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class HashBagTest extends Specification {

    def 'should be empty on construction'() {
        expect:
        new HashBag().isEmpty()
        new HashBag(16).isEmpty()
        new HashBag([]).isEmpty()
    }

    def 'should add elements and track size'() {
        given:
        def bag = new HashBag()

        when:
        bag.add('a')
        bag.add('b')
        bag.add('a')

        then:
        bag.size() == 3
        !bag.isEmpty()
    }

    def 'should construct from collection'() {
        given:
        def bag = new HashBag(['x', 'y', 'x', 'z'])

        expect:
        bag.size() == 4
        bag.contains('x')
        bag.contains('y')
        bag.contains('z')
    }

    def 'should contain added element'() {
        given:
        def bag = new HashBag()

        when:
        bag.add('hello')

        then:
        bag.contains('hello')
        !bag.contains('world')
    }

    def 'should track duplicate counts'() {
        given:
        def bag = new HashBag()

        when:
        bag.add(1)
        bag.add(1)
        bag.add(1)

        then:
        bag.size() == 3
        bag.contains(1)
    }

    def 'should addAll elements from a collection'() {
        given:
        def bag = new HashBag()

        when:
        bag.addAll([1, 2, 2, 3])

        then:
        bag.size() == 4
        bag.contains(1)
        bag.contains(2)
        bag.contains(3)
    }

    def 'should remove a single instance of an element'() {
        given:
        def bag = new HashBag(['a', 'a', 'b'])

        when:
        def removed = bag.remove('a')

        then:
        removed
        bag.size() == 2
        bag.contains('a')
        bag.contains('b')
    }

    def 'should remove last instance and stop containing element'() {
        given:
        def bag = new HashBag(['a', 'b'])

        when:
        bag.remove('a')

        then:
        !bag.contains('a')
        bag.contains('b')
        bag.size() == 1
    }

    def 'should return false when removing absent element'() {
        given:
        def bag = new HashBag(['a'])

        expect:
        !bag.remove('z')
        bag.size() == 1
    }

    def 'should removeAll one instance of each listed element'() {
        given:
        def bag = new HashBag(['a', 'a', 'b', 'c'])

        when:
        def changed = bag.removeAll(['a', 'c'])

        then:
        changed
        bag.size() == 2
        bag.contains('a')
        !bag.contains('c')
    }

    def 'should clear all elements'() {
        given:
        def bag = new HashBag([1, 2, 3])

        when:
        bag.clear()

        then:
        bag.isEmpty()
        bag.size() == 0
    }

    def 'should containsAll when all elements present'() {
        given:
        def bag = new HashBag(['a', 'b', 'c'])

        expect:
        bag.containsAll(['a', 'b'])
        bag.containsAll(['a', 'b', 'c'])
        !bag.containsAll(['a', 'b', 'c', 'd'])
    }

    def 'should retainAll only elements in the given collection'() {
        given:
        def bag = new HashBag(['a', 'b', 'c', 'a'])

        when:
        def changed = bag.retainAll(['a', 'c'])

        then:
        changed
        bag.contains('a')
        bag.contains('c')
        !bag.contains('b')
    }

    def 'should return false from retainAll when nothing removed'() {
        given:
        def bag = new HashBag(['a', 'b'])

        when:
        def changed = bag.retainAll(['a', 'b', 'c'])

        then:
        !changed
        bag.size() == 2
    }

    def 'should iterate over all elements including duplicates'() {
        given:
        def bag = new HashBag(['x', 'x', 'y'])

        when:
        def result = []
        for( def e : bag )
            result << e

        then:
        result.size() == 3
        result.count('x') == 2
        result.count('y') == 1
    }

    def 'should iterate over empty bag without error'() {
        given:
        def bag = new HashBag()

        expect:
        !bag.iterator().hasNext()
    }

    def 'should be equal to bag with same element counts'() {
        given:
        def bag1 = new HashBag(['a', 'b', 'a'])
        def bag2 = new HashBag(['a', 'b', 'a'])
        def bag3 = new HashBag(['a', 'b'])

        expect:
        bag1 == bag2
        bag1 != bag3
    }

    def 'should not be equal to non-HashBag object'() {
        given:
        def bag = new HashBag(['a'])

        expect:
        bag != ['a']
        bag != null
    }

    def 'should have consistent hashCode with equals'() {
        given:
        def bag1 = new HashBag(['a', 'b', 'a'])
        def bag2 = new HashBag(['a', 'b', 'a'])

        expect:
        bag1.hashCode() == bag2.hashCode()
    }

}
