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

package nextflow.splitter

import nextflow.Channel
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BytesSplitterTest extends Specification {

    def bytes = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6 ] as byte[]

    def testSplitterCount() {
        expect:
        new BytesSplitter().options(by: 6).target(bytes).count() == 3

    }


    def testSplitterList() {

        expect:
        new BytesSplitter().options(by: 5).target(bytes).list() == [ [0, 1, 2, 3, 4] as byte[], [5, 6, 7, 8, 9] as byte[] , [ 0, 1, 2, 3, 4] as byte[], [5, 6] as byte[] ]

    }

    def testSplitterWithLimit() {

        expect:
        new BytesSplitter().options(by: 5, limit: 10).target(bytes).list() == [ [0, 1, 2, 3, 4] as byte[], [5, 6, 7, 8, 9] as byte[] ]
        new BytesSplitter().options(by: 5, limit: 13).target(bytes).list() == [ [0, 1, 2, 3, 4] as byte[], [5, 6, 7, 8, 9] as byte[], [ 0, 1, 2 ] as byte[] ]

    }

    def testSplitterChannel() {

        when:
        def c = new BytesSplitter().options(by: 5).target(bytes).channel()
        then:
        c.val == [0, 1, 2, 3, 4] as byte[]
        c.val == [5, 6, 7, 8, 9] as byte[]
        c.val == [ 0, 1, 2, 3, 4] as byte[]
        c.val == [5, 6] as byte[]
        c.val == Channel.STOP

    }

}
