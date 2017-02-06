/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractSplitterTest extends Specification {


    def 'test invokeEachClosure method'() {

        given:
        def file = Paths.get('/some/file.txt')

        when:
        def splitter = [:] as AbstractSplitter
        then:
        splitter.invokeEachClosure(null, 'hola') == 'hola'
        splitter.invokeEachClosure({ x -> x.reverse() }, 'hola') == 'aloh'

        when:
        splitter.elem = 0
        splitter.target( [1,2,3] )
        then:
        splitter.invokeEachClosure(null, 'hola') == ['hola', 2, 3]
        splitter.invokeEachClosure( { x -> x } , 'hola') == ['hola', 2, 3]
        splitter.invokeEachClosure( { a,b,c -> a } , 'hola') == 'hola'
        splitter.invokeEachClosure( { tuple -> tuple[0] } , 'hola') == 'hola'
    }

    def 'test isTrueOrMap'() {
        expect:
        AbstractSplitter.isTrueOrMap(true)
        AbstractSplitter.isTrueOrMap(Boolean.TRUE)
        AbstractSplitter.isTrueOrMap([:])
        !AbstractSplitter.isTrueOrMap(0)
        !AbstractSplitter.isTrueOrMap(false)
        !AbstractSplitter.isTrueOrMap(Boolean.FALSE)

    }

    def 'test findSource'() {

        def splitter

        when:
        splitter = [:] as AbstractSplitter
        then:
        splitter.elem == -1
        splitter.findSource([ 10, 20 ]) == 10
        splitter.elem == 0

        when:
        splitter = [:] as AbstractSplitter
        then:
        splitter.elem == -1
        splitter.findSource([ 10, Paths.get('/hello') ]) == Paths.get('/hello')
        splitter.elem == 1


        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = 1
        then:
        splitter.findSource([ 10, 20, Paths.get('/hello') ]) == 20


    }

}
