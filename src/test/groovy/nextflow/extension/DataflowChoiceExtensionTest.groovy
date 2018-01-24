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
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowChoiceExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }


    @Timeout(1)
    def 'should choice targets specified with an open array'() {

        given:
        def source = Channel.from 'Hello world', 'Hola', 'Hello John'
        def queue1 = Channel.create()
        def queue2 = Channel.create()
        def queue3 = Channel.create()

        when:
        source.choice( queue1, queue2, queue3 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue1.val == 'Hello John'
        queue1.val == Channel.STOP
        queue2.val == 'Hola'
        queue2.val == Channel.STOP
        queue3.val == Channel.STOP

    }

    def 'should choice targets specified with a list '() {

        given:
        def source = Channel.from 'Hello world', 'Hola', 'Hello John'
        def queue1 = Channel.create()
        def queue2 = Channel.create()

        when:
        source.choice( [queue1, queue2] ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue1.val == 'Hello John'
        queue1.val == Channel.STOP
        queue2.val == 'Hola'
        queue2.val == Channel.STOP

    }


    @Timeout(1)
    def 'should choice from a dataflow value channel'() {

        given:
        def source = Channel.value('Hello world')
        def queue1 = new DataflowVariable()
        def queue2 = new DataflowVariable()

        when:
        source.choice( queue1, queue2 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

        then:
        queue1.val == 'Hello world'
        queue2.val == Channel.STOP

    }


}
