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

import nextflow.Global
import nextflow.Session
import nextflow.executor.Executor
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskRunTest2 extends Specification {

    @Unroll
    def 'should get container platform' () {
        given:
        def session = Global.session = Mock(Session) { getConfig()>>SESSION }
        def executor = Mock(Executor) { getSession()>>session }
        def processor = Mock(TaskProcessor) { getExecutor()>>executor; getSession()>>session }
        and:
        def config = new TaskConfig(CONFIG)
        def task = new TaskRun(config: config, processor: processor)

        expect:
        task.getContainerPlatform() == EXPECTED

        where:
        CONFIG          | SESSION                   | EXPECTED
        [:]             | [:]                       | null
        [arch:'amd64']  | [:]                       | 'linux/amd64'
        [arch:'arm64']  | [:]                       | 'linux/arm64'
        and:
        [:]             | [wave:[enabled:true]]     | 'linux/amd64'
        [arch:'amd64']  | [wave:[enabled:true]]     | 'linux/amd64'
        [arch:'arm64']  | [wave:[enabled:true]]     | 'linux/arm64'
    }

}
