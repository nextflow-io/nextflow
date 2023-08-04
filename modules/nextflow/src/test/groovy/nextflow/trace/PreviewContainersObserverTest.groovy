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

package nextflow.trace

import nextflow.dag.DAG
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PreviewContainersObserverTest extends Specification {

    def makeVertex(DAG dag, String name, String container) {
        final processor = Mock(TaskProcessor)
        processor.name >> name
        processor.getPreviewTask() >> Mock(TaskRun) {
            getContainer() >> container
        }

        final vertex = new DAG.Vertex(dag, DAG.Type.PROCESS)
        vertex.process = processor

        return vertex
    }

    def 'should get containers' () {
        given:
        def dag = Mock(DAG)
        dag.vertices >> [
            makeVertex(dag, 'proc1', 'container1'),
            makeVertex(dag, 'proc2', 'container2')
        ]

        when:
        def observer = new PreviewContainersObserver(dag: dag)
        then:
        observer.getContainers() == [
            'proc1': 'container1',
            'proc2': 'container2',
        ]
    }

    def 'should render config output' () {
        when:
        def observer = new PreviewContainersObserver()
        def containers = [
            'proc1': 'container1',
            'proc2': 'container2',
        ]
        then:
        observer.renderConfig(containers) == '''\
            process { withName: 'proc1' { container = 'container1' } }
            process { withName: 'proc2' { container = 'container2' } }
            '''.stripIndent()
    }

    def 'should render json output' () {
        when:
        def observer = new PreviewContainersObserver()
        def containers = [
            'proc1': 'container1',
            'proc2': 'container2',
        ]
        then:
        observer.renderJson(containers) == '''\
            [
                {
                    "name": "proc1",
                    "container": "container1"
                },
                {
                    "name": "proc2",
                    "container": "container2"
                }
            ]
            '''.stripIndent().trim()
    }

}
