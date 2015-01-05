/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import nextflow.script.TaskBody
import nextflow.util.BlankSeparatedList
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskContextTest extends Specification {

    def testSaveAndReadContextMap () {

        setup:
        def taskConfig = new ProcessConfig([:])
        def file = Files.createTempFile('test.ctx',null)
        def processor = [:] as TaskProcessor
        processor.metaClass.getTaskConfig = { taskConfig }
        processor.metaClass.getTaskBody = { new TaskBody(null,'source',true) }
        def str = 'Hola'
        def map = new TaskContext(processor, [:])
        map.alpha = 1
        map.beta = "${str}.txt"
        map.delta = new Duration('1day')
        map.micro = new MemoryUnit('100KB')
        map.file = Paths.get('Hola.txt')
        map.list = new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        map.holder = 'just a string'
        map.context = [uno: 1, due: 'str', tre: "$str"]

        when:
        map.save(file)
        def result = TaskContext.read(processor, file)

        then:
        result.size() == 8
        result.alpha == 1
        result.beta == "${str}.txt"
        result.delta == new Duration('1day')
        result.micro == new MemoryUnit('100KB')
        result.file.equals( Paths.get('Hola.txt') )
        result.context == [uno: 1, due: 'str', tre: "$str"]
        result.list == new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        result.holder == 'just a string'
        result.get('holder') == 'just a string'
        result.getHolder() instanceof Map
        result.getHolder().get('alpha') == 1

        cleanup:
        file.delete()

    }


    def testDehydrateRehydrate() {

        setup:
        def bind = new Binding(x:1, y:2)
        def script = new Script() {
            @Override
            Object run() {
                return null
            }
        }
        script.setBinding(bind)

        def local = [p:3, q:4, path: Paths.get('some/path')]
        def delegate = new TaskContext( script, local, 'hola' )

        when:
        def bytes = delegate.dehydrate()
        def copy = TaskContext.rehydrate(bytes)

        then:
        delegate == copy
        delegate.getHolder() == copy.getHolder()
        copy.getHolder() == local

    }

}
