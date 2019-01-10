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
