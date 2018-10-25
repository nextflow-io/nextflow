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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BlankSeparatedListTest extends Specification {

    def 'should return a blank separated list of values'() {

        expect:
        new BlankSeparatedList('a'..'z').toString() == ('a'..'z').join(' ')
        "${new BlankSeparatedList('a'..'z')}" == ('a'..'z').join(' ')

    }


    def 'should access entry with square brackets'() {

        when:
        def x = new BlankSeparatedList('a'..'z')
        then:
        x[0] == 'a'
        x[1] == 'b'
        x[2] == 'c'
    }

    def 'should collect items' () {

        given:
        def p = new BlankSeparatedList('x'..'z')
        expect:
        p.collect() == ['x','y','z']

    }

    def 'should convert to a list' () {
        given:
        def p = new BlankSeparatedList('x'..'z')
        expect:
        p as List == ['x','y','z']
        p as Set == ['x','y','z'] as Set
    }

    def 'should join items' () {
        given:
        def p = new BlankSeparatedList('x'..'z')
        expect:
        p.join('-') == 'x-y-z'
    }

    def 'should get first item' () {
        given:
        def p = new BlankSeparatedList('x'..'z')
        expect:
        p.first() == 'x'
    }


    def 'should serialised blank separated list' () {

        given:
        def list = new BlankSeparatedList(['alpha','beta',null,'gamma'])
        when:
        def buffer = KryoHelper.serialize(list)
        then:
        KryoHelper.deserialize(buffer) == list
        KryoHelper.deserialize(buffer) instanceof BlankSeparatedList

    }

}
