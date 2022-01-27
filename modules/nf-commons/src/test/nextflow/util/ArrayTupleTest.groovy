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
class ArrayTupleTest extends Specification {

    def testTuple() {
        when:
        def tuple = new ArrayTuple([1,2,3])
        then:
        tuple[0] == 1
        tuple[1] == 2
        tuple[2] == 3
        tuple.size() == 3

        when:
        tuple.add(4)
        then:
        thrown(UnsupportedOperationException)

        when:
        tuple.remove(0)
        then:
        thrown(UnsupportedOperationException)

        when:
        tuple[0] = 2
        then:
        thrown(UnsupportedOperationException)

        expect:
        tuple.sum() == 6
        tuple.iterator().next() == 1


    }
}
