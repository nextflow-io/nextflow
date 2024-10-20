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
class TaskDAGTest extends Specification {

    def 'should add task vertices and outputs' () {

        given:
        def task1 = Mock(TaskRun) {
            getInputFilesMap() >> [
                'data.txt': Paths.get('/inputs/data.txt')
            ]
            getOutputsByType(_) >> [
                'data.foo': Paths.get('/work/00112233/data.foo')
            ]
        }
        def task2 = Mock(TaskRun) {
            getInputFilesMap() >> [
                'data.foo': Paths.get('/work/00112233/data.foo')
            ]
            getOutputsByType(_) >> [
                'data.bar': Paths.get('/work/aabbccdd/data.bar')
            ]
        }
        def dag = new TaskDAG()

        when:
        dag.addTask( task1 )
        dag.addTask( task2 )
        def v1 = dag.vertices[task1]
        def v2 = dag.vertices[task2]
        then:
        v1.inputs.size() == 1
        v1.inputs['data.txt'] == Paths.get('/inputs/data.txt')
        and:
        v2.inputs.size() == 1
        v2.inputs['data.foo'] == Paths.get('/work/00112233/data.foo')

        when:
        dag.addTaskOutputs( task1 )
        dag.addTaskOutputs( task2 )
        then:
        v1.outputs == [ Paths.get('/work/00112233/data.foo') ] as Set
        and:
        v2.outputs == [ Paths.get('/work/aabbccdd/data.bar') ] as Set
        and:
        dag.getProducerVertex(v1.inputs['data.txt']) == null
        dag.getProducerVertex(v2.inputs['data.foo']) == v1
    }

}
