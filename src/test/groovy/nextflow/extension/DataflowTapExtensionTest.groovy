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
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowTapExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }


    def 'should `tap` item to a new channel' () {

        when:
        def result = Channel.from( 4,7,9 ) .tap { first }.map { it+1 }
        then:
        session.binding.first.val == 4
        session.binding.first.val == 7
        session.binding.first.val == 9
        session.binding.first.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` item to more than one channel' () {

        when:
        def result = Channel.from( 4,7,9 ) .tap { foo; bar }.map { it+1 }
        then:
        session.binding.foo.val == 4
        session.binding.foo.val == 7
        session.binding.foo.val == 9
        session.binding.foo.val == Channel.STOP
        session.binding.bar.val == 4
        session.binding.bar.val == 7
        session.binding.bar.val == 9
        session.binding.bar.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` target channel' () {

        when:
        def target = Channel.create()
        def result = Channel.from( 8,2,5 ) .tap(target).map { it+1 }
        then:
        result instanceof DataflowQueue
        target instanceof DataflowQueue

        target.val == 8
        target.val == 2
        target.val == 5
        target.val == Channel.STOP

        result.val == 9
        result.val == 3
        result.val == 6
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` dataflow value' () {

        when:
        def target = Channel.value()
        def result = Channel.value(7) .tap(target).map { it+1 }
        then:
        result instanceof DataflowVariable
        target instanceof DataflowVariable

        target.val == 7
        target.val == 7
        result.val == 8
        result.val == 8

        !session.dag.isEmpty()

    }

    def 'should `tap` dataflow value and target as queue' () {

        when:
        new DataflowVariable() .tap( new DataflowQueue() )
        then:
        thrown(IllegalArgumentException)
        session.dag.isEmpty()

    }

}
