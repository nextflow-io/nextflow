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

package nextflow.processor

import nextflow.util.KryoHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskIdTest extends Specification {

    def 'should validate equals' () {

        expect:
        new TaskId(1) == new TaskId(1)
        new TaskId(1) != new TaskId(2)
        new TaskId(1) == 1
        new TaskId(1) != 2
        1 == new TaskId(1)
        1 != new TaskId(2)

    }

    def 'should validate comparator' () {

        expect:
        new TaskId(1) < new TaskId(9)
        new TaskId(1) <= new TaskId(9)
        new TaskId(5) > new TaskId(2)
        new TaskId(5) >= new TaskId(2)

        new TaskId(1) < 9
        new TaskId(1) <= 9
        new TaskId(5) > 2
        new TaskId(5) >= 2

        1 < new TaskId(9)
        1 <= new TaskId(9)
        5 > new TaskId(2)
        5 >= new TaskId(2)

        new TaskId(0) < Integer.MAX_VALUE
        Integer.MIN_VALUE < new TaskId(0)

    }

    def 'should validate toString' () {

        expect:
        new TaskId(1).toString() == '1'
        new TaskId(100).toString() == '100'

    }

    def 'should validate factory method' () {
        expect:
        TaskId.of(100) == new TaskId(100)
        TaskId.of('100') == new TaskId(100)
        TaskId.of('101') != new TaskId(100)
    }

    def 'should serialise/deserialise key objects' () {

        given:
        final keys = [ TaskId.of(10), TaskId.of(20) ]

        when:
        def bytes = KryoHelper.serialize(keys)
        then:
        ((List)KryoHelper.deserialize(bytes)) == keys
    }

}
