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

package nextflow.extension
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowIntoExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }


    def 'should create two target channels'() {
        given:
        def result = Channel.from(1,2,3,4)

        when:
        def (ch1, ch2) = result.into(2)

        then:
        ch1 instanceof DataflowQueue
        ch2 instanceof DataflowQueue

        ch1.val == 1
        ch1.val == 2
        ch1.val == 3
        ch1.val == 4
        ch1.val == Channel.STOP

        ch2.val == 1
        ch2.val == 2
        ch2.val == 3
        ch2.val == 4
        ch2.val == Channel.STOP

    }

    def 'should forward items into target channels'() {
        given:
        def result = Channel.from('a','b',[1,2])
        def ch1 = Channel.create()
        def ch2 = Channel.create()
        def ch3 = Channel.create()

        when:
        result.into(ch1, ch2, ch3)

        then:
        ch1.val == 'a'
        ch1.val == 'b'
        ch1.val == [1,2]
        ch1.val == Channel.STOP

        ch2.val == 'a'
        ch2.val == 'b'
        ch2.val == [1,2]
        ch2.val == Channel.STOP

        ch3.val == 'a'
        ch3.val == 'b'
        ch3.val == [1,2]
        ch3.val == Channel.STOP

    }

    @Timeout(1)
    def 'should forward dataflow value into a new channel'() {

        when:
        def result = Channel.create()
        Channel.value('Hello').into(result)
        then:
        result.val == 'Hello'
        result.val == Channel.STOP

    }

    def 'should create new dataflow variables and forward item to them'  () {

        when:
        Channel.from(10,2,30).into { alpha; gamma }

        then:
        session.binding.alpha.val == 10
        session.binding.alpha.val == 2
        session.binding.alpha.val == 30
        session.binding.alpha.val == Channel.STOP

        session.binding.gamma.val == 10
        session.binding.gamma.val == 2
        session.binding.gamma.val == 30
        session.binding.gamma.val == Channel.STOP

    }


}
