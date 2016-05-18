/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EachInParamTest extends Specification {

    def testNormalize() {

        given:
        def channel = Channel.from(1,2,3,5)
        def value = Channel.value('a')
        def list = Channel.value([4,5,6])

        expect:
        EachInParam.normalizeToVariable(1).getVal() == [1]
        EachInParam.normalizeToVariable([3,4,5]).getVal() == [3,4,5]
        EachInParam.normalizeToVariable(channel).get() == [1,2,3,5]
        EachInParam.normalizeToVariable(value).get() == ['a']
        EachInParam.normalizeToVariable(list).get() == [4,5,6]





    }

}
