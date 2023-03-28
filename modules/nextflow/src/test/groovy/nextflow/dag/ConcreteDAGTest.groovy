/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.dag

import java.nio.file.Paths

import com.google.common.hash.HashCode
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConcreteDAGTest extends Specification {

    def 'should add task nodes and outputs' () {

        given:
        def task1 = Mock(TaskRun) {
            getHash() >> HashCode.fromString('00112233445566778899aabbccddeeff')
            getName() >> 'foo'
            getInputFilesMap() >> [
                'data.txt': Paths.get('/inputs/data.txt')
            ]
            getOutputsByType(_) >> [
                'data.foo': Paths.get('/work/00/112233445566778899aabbccddeeff/data.foo')
            ]
        }
        def task2 = Mock(TaskRun) {
            getHash() >> HashCode.fromString('aabbccddeeff00112233445566778899')
            getName() >> 'bar'
            getInputFilesMap() >> [
                'data.foo': Paths.get('/work/00/112233445566778899aabbccddeeff/data.foo')
            ]
            getOutputsByType(_) >> [
                'data.bar': Paths.get('/work/aa/bbccddeeff00112233445566778899/data.bar')
            ]
        }
        def dag = new ConcreteDAG()

        when:
        dag.addTask( task1 )
        dag.addTask( task2 )
        def node1 = dag.nodes['00112233445566778899aabbccddeeff']
        def node2 = dag.nodes['aabbccddeeff00112233445566778899']
        then:
        node1.index == 0
        node1.label == '[00/112233] foo'
        node1.inputs.size() == 1
        node1.inputs[0].name == 'data.txt'
        node1.inputs[0].path == Paths.get('/inputs/data.txt')
        node1.inputs[0].predecessor == null
        node2.index == 1
        node2.label == '[aa/bbccdd] bar'
        node2.inputs.size() == 1
        node2.inputs[0].name == 'data.foo'
        node2.inputs[0].path == Paths.get('/work/00/112233445566778899aabbccddeeff/data.foo')
        node2.inputs[0].predecessor == '00112233445566778899aabbccddeeff'

        when:
        dag.addTaskOutputs( task1 )
        dag.addTaskOutputs( task2 )
        then:
        node1.outputs.size() == 1
        node1.outputs[0].name == 'data.foo'
        node1.outputs[0].path == Paths.get('/work/00/112233445566778899aabbccddeeff/data.foo')
        node2.outputs.size() == 1
        node2.outputs[0].name == 'data.bar'
        node2.outputs[0].path == Paths.get('/work/aa/bbccddeeff00112233445566778899/data.bar')
    }

}
