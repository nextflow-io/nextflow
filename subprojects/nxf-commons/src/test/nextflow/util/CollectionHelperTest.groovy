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
class CollectionHelperTest extends Specification {

    def testFlatten() {

        expect:
        CollectionHelper.flatten( [1,2,3] )  == [ [1,2,3] ]

        CollectionHelper.flatten( ['a', [1,2,3]] ) == [ ['a', 1], ['a', 2], ['a', 3] ]

        CollectionHelper.flatten( [ ['a','b'], [1,2,3,4]] ) == [
                ['a', 1],
                ['b', 2],
                [null, 3],
                [null, 4],

        ]


        CollectionHelper.flatten( [ 'x',  ['a','b'], [1,2,3,4]] ) == [
                ['x', 'a', 1],
                ['x', 'b', 2],
                ['x', null, 3],
                ['x', null, 4],
        ]

    }


}
