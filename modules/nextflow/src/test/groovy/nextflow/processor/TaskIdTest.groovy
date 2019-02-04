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
