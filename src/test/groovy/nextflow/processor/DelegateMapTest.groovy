/*
 * Copyright (c) 2012, the authors.
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
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.script.BaseScript
import nextflow.util.BlankSeparatedList
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DelegateMapTest extends Specification {

    def 'test undef'( ) {

        setup:
        def script = Mock(BaseScript)
        def map = new DelegateMap(script)

        when:
        map.x = 1
        then:
        map.x == 1

        when:
        map.get('y')
        then:
        thrown(MissingPropertyException)

        when:
        def val = new DelegateMap(script,true).get('y')
        then:
        val == '$y'


    }


    def testSaveAndReadContextMap () {

        setup:
        def taskConfig = new TaskConfig([undef:false])
        def file = Files.createTempFile('test.ctx',null)
        def processor = [:] as TaskProcessor
        processor.metaClass.getTaskConfig = { taskConfig }
        def str = 'Hola'
        def map = new DelegateMap(processor)
        map.alpha = 1
        map.beta = "${str}.txt"
        map.file = Paths.get('Hola.txt')
        map.list = new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        map.holder = 'just a string'

        when:
        map.save(file)
        def result = DelegateMap.read(processor, file)

        then:
        result.size() == 5
        result.alpha == 1
        result.beta == "${str}.txt"
        result.file.equals( Paths.get('Hola.txt') )
        result.list == new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        result.holder == 'just a string'
        result.get('holder') == 'just a string'
        result.getHolder() instanceof Map
        result.getHolder().get('alpha') == 1

        cleanup:
        file.delete()

    }




}
