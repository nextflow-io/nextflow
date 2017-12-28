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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ObjectListCollectorTest extends Specification {

    def 'test list text buffer' () {

        given:
        def buffer = new ObjectListCollector()

        when:
        buffer.add('alpha')
        buffer.add('delta')
        buffer.add('gamma')

        then:
        buffer.nextChunk() == ['alpha','delta','gamma']

    }


    def 'test is empty' () {

        when:
        def buffer = new ObjectListCollector()
        then:
        !buffer.hasChunk()

        when:
        buffer.add('hello')
        then:
        buffer.hasChunk()

        when:
        buffer.nextChunk()
        then:
        !buffer.hasChunk()

    }

    def 'test reset' () {

        given:
        def buffer = new ObjectListCollector()
        buffer.add('hello')

        when:
        buffer.nextChunk()
        then:
        !buffer.hasChunk()
        buffer.nextChunk() == []

    }

}
