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

    def 'should add task vertices and outputs' () {

        given:
        def task1 = Mock(TaskRun) {
            getHash() >> HashCode.fromString('00112233')
            getName() >> 'foo'
            getInputFilesMap() >> [
                'data.txt': Paths.get('/inputs/data.txt')
            ]
            getOutputsByType(_) >> [
                'data.foo': Paths.get('/work/00112233/data.foo')
            ]
        }
        def task2 = Mock(TaskRun) {
            getHash() >> HashCode.fromString('aabbccdd')
            getName() >> 'bar'
            getInputFilesMap() >> [
                'data.foo': Paths.get('/work/00112233/data.foo')
            ]
            getOutputsByType(_) >> [
                'data.bar': Paths.get('/work/aabbccdd/data.bar')
            ]
        }
        def dag = new ConcreteDAG()

        when:
        dag.addTask( task1 )
        dag.addTask( task2 )
        def v1 = dag.vertices[task1]
        def v2 = dag.vertices[task2]
        then:
        v1.index == 0
        v1.label == '[00/112233] foo'
        v1.inputs.size() == 1
        v1.inputs[0] == Paths.get('/inputs/data.txt')
        and:
        v2.index == 1
        v2.label == '[aa/bbccdd] bar'
        v2.inputs.size() == 1
        v2.inputs[0] == Paths.get('/work/00112233/data.foo')

        when:
        dag.addTaskOutputs( task1 )
        dag.addTaskOutputs( task2 )
        then:
        v1.outputs.size() == 1
        v1.outputs[0] == Paths.get('/work/00112233/data.foo')
        and:
        v2.outputs.size() == 1
        v2.outputs[0] == Paths.get('/work/aabbccdd/data.bar')
        and:
        dag.getProducerVertex(v1.inputs[0]) == null
        dag.getProducerVertex(v2.inputs[0]) == v1
    }

    def 'should write meta file' () {

        given:
        def folder = File.createTempDir()
        def outputFile = new File(folder, 'data.bar') ; outputFile.text = 'bar'

        def task1 = Mock(TaskRun) {
            name >> 'foo'
            hash >> HashCode.fromString('00112233')
            getInputFilesMap() >> [ 'data.txt': Paths.get('/inputs/data.txt') ]
            getOutputsByType(_) >> [ 'data.foo': Paths.get('/work/00112233/data.foo') ]
        }
        def task2 = Mock(TaskRun) {
            name >> 'bar'
            workDir >> folder.toPath()
            hash >> HashCode.fromString('aabbccdd')
            getInputFilesMap() >> [ 'data.foo': Paths.get('/work/00112233/data.foo') ]
            getOutputsByType(_) >> [ 'data.bar': outputFile.toPath() ]
        }
        def dag = new ConcreteDAG()

        when:
        dag.addTask(task1)
        dag.addTaskOutputs(task1)
        dag.addTask(task2)
        dag.addTaskOutputs(task2)
        dag.writeMetaFile(task2)
        then:
        task2.workDir.resolve(TaskRun.CMD_META).text == """{"hash":"aabbccdd","inputs":[{"name":"data.foo","path":"/work/00112233/data.foo","predecessor":"00112233"}],"outputs":[{"name":"data.bar","path":"${folder}/data.bar","size":3,"checksum":"37b51d194a7513e45b56f6524f2d51f2"}]}"""

        cleanup:
        folder.delete()
    }

}
