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

import java.nio.file.Path
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractSplitterTest extends Specification {

    private Path file(String name) { Paths.get(name) }

    def 'test invokeEachClosure method'() {


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
        splitter.elem == null
        splitter.findSource([ 10, 20 ]) == 10
        splitter.elem == 0

        when:
        splitter = [:] as AbstractSplitter
        then:
        splitter.elem == null
        splitter.findSource([ 10, file('/hello') ]) == file('/hello')
        splitter.elem == 1

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = -2  // <-- find the second file
        then:
        splitter.findSource([ 10, 20, file('/hello'), 30, 40, file('/second'), file('/third') ]) == file('/second')
        splitter.elem == 5

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = -3  // <-- find the third file
        then:
        splitter.findSource([ 10, 20, file('/hello'), 30, 40, file('/second'), file('/third') ]) == file('/third')
        splitter.elem == 6

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = 1
        then:
        splitter.findSource([ 10, 20, file('/hello') ]) == 20
        splitter.elem == 1

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = 0
        then:
        splitter.findSource([ 10, 20, Paths.get('/hello') ]) == 10
        splitter.elem == 0

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = -1
        splitter.findSource([ 10, 20 ])
        then:
        thrown(IllegalArgumentException)

        when:
        splitter = [:] as AbstractSplitter
        splitter.elem = -2
        splitter.findSource([ 10, 20, Paths.get('/hello') ])
        then:
        thrown(IllegalArgumentException)
    }

}
