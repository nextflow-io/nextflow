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

package nextflow.splitter

import groovyx.gpars.dataflow.operator.PoisonPill
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StringSplitterTest extends Specification {


    def testSplitString () {

        expect:
        new StringSplitter().options(by: 5).target('012345678901234567') .list() == ['01234','56789','01234','567']

    }

    def testSplitStringWithLimit () {

        expect:
        new StringSplitter().options(by:5, limit: 11).target('012345678901234567') .list() == ['01234','56789','0']

    }

    def testSplitWithClosure() {

        expect:
        new StringSplitter()
            .target('012345678901234567')
            .options(by:5, each: {it.reverse()} )
            .list()  == ['43210','98765','43210','765']

    }

    def testSplitChannel() {

        when:
        def q = new StringSplitter().target('012345678901234567') .options(by:5). channel()
        then:
        q.val == '01234'
        q.val == '56789'
        q.val == '01234'
        q.val == '567'
        q.val == PoisonPill.instance

    }

    def testSplitStringByOne () {

        expect:
        new StringSplitter().options(by: 1).target('ABC') .list() == ['A','B','C']

    }

    def testSplitToFile() {
        given:
        def folder = TestHelper.createInMemTempDir()

        when:
        def chunks = new StringSplitter().options(by:4, file: folder).target('Hello world') .list()
        then:
        chunks[0].text == 'Hell'
        chunks[1].text == 'o wo'
        chunks[2].text == 'rld'

    }

}
