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

package nextflow.script

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EachInParamTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def testNormalize() {

        given:
        def channel = Channel.from(1,2,3,5)
        def value = Channel.value('a')
        def list = Channel.value([4,5,6])
        def each = new EachInParam(Mock(Binding), [])

        expect:
        each.normalizeToVariable(1).val == [1]
        each.normalizeToVariable([3,4,5]).val == [3,4,5]
        each.normalizeToVariable(channel).val == [1,2,3,5]
        each.normalizeToVariable(value).val == ['a']
        each.normalizeToVariable(list).val == [4,5,6]

    }

}
