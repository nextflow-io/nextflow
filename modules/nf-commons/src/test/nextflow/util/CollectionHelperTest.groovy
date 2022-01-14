/*
 * Copyright 2020-2022, Seqera Labs
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
