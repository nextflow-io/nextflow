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
