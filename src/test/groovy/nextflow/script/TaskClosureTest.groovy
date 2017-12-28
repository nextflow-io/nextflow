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

package nextflow.script

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskClosureTest extends Specification {

    def 'should invoke inner closure' () {

        given:
        def code = { String it -> it.reverse() }

        when:
        def wrap = new TaskClosure(code, 'reverse')
        then:
        wrap.getMaximumNumberOfParameters() == 1
        wrap.getParameterTypes() as List == [String]
        wrap.call('Hello') == 'olleH'

    }

    def 'should clone and invoke the closure' ()  {

        given:
        def code = { String a, Integer b -> a.toInteger() + b }
        def wrap = new TaskClosure(code, 'closure source code')

        when:
        def copy = (TaskClosure)wrap.clone()
        then:
        copy.getSource() == 'closure source code'
        !copy.getInnerClosure().is( code )
        copy.getMaximumNumberOfParameters() == 2
        copy.getParameterTypes() as List == [String, Integer]
        copy.call('3', 2) == 5

    }
}
