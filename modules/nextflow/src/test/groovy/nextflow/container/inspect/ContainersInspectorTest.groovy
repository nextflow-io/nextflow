/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.container.inspect

import nextflow.NF
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.ProcessDef
import nextflow.script.ScriptMeta
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ContainersInspectorTest extends Specification {

    def setup() {
        ScriptMeta.reset()
        NF.init()
    }

    def makeProcess(String name, String container) {
        final script = new BaseScript() {
            Object runScript() {}
        }
        final processDef = Mock(ProcessDef) {
            getName() >> name
            createTaskProcessor() >> Mock(TaskProcessor) {
                createTaskPreview() >> Mock(TaskRun) {
                    getContainer() >> container
                }
            }
        }
        ScriptMeta.get(script).addDefinition(processDef)
    }

    def 'should get containers' () {
        given:
        makeProcess('proc1', 'container1')
        makeProcess('proc2', 'container2')

        when:
        def observer = new ContainersInspector(false)
        then:
        observer.getContainers() == [
            'proc1': 'container1',
            'proc2': 'container2',
        ]
    }

    def 'should render containers as json' () {
        given:
        makeProcess('proc1', 'container1')
        makeProcess('proc2', 'container2')

        when:
        def result = new ContainersInspector(false)
                .withFormat('json')
                .renderContainers()
        then:
        result == '''\
            {
                "processes": [
                    {
                        "name": "proc1",
                        "container": "container1"
                    },
                    {
                        "name": "proc2",
                        "container": "container2"
                    }
                ]
            }
            '''.stripIndent(true)
    }

    def 'should render containers as nextflow config' () {
        given:
        makeProcess('proc1', 'container1')
        makeProcess('proc2', 'container2')

        when:
        def result = new ContainersInspector(false)
                .withFormat('config')
                .renderContainers()
        then:
        result == '''\
            process { withName: 'proc1' { container = 'container1' } }
            process { withName: 'proc2' { container = 'container2' } }
            '''.stripIndent(true)
    }

}
