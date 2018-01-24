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

}
