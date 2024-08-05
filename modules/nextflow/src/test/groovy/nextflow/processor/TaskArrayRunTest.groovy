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

package nextflow.processor

import nextflow.Session
import nextflow.container.resolver.ContainerInfo
import nextflow.executor.Executor
import nextflow.executor.TaskArrayExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskArrayRunTest extends Specification {

    static abstract class TestExecutor extends Executor implements TaskArrayExecutor {

    }

    def 'should get container info' () {
        given:
        def session = Mock(Session)
        def executor = Mock(TestExecutor) { getSession()>>session }
        def processor = Mock(TaskProcessor) { getExecutor()>>executor; getSession()>>session }
        and:
        def config = new TaskConfig([container:'ubuntu'])
        def task = new TaskArrayRun(config: config, processor: processor)
        when:
        def info = task.containerInfo()
        then:
        info == new ContainerInfo('ubuntu','ubuntu','ubuntu')
    }

    def 'should be an array' () {
        given:
        def task = new TaskArrayRun()
        expect:
        task.isArray()
    }

}
