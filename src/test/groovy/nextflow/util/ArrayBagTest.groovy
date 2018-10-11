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

}
