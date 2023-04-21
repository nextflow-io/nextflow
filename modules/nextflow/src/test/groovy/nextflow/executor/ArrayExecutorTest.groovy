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

package nextflow.executor

import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ArrayExecutorTest extends Specification {

    def 'should throw error if target executor does not support array jobs' () {

        setup:
        def executorFactory = Spy(ExecutorFactory)
        def session = Mock(Session) {
            getExecutorFactory() >> executorFactory
        }
        def executor = new ArrayExecutor()
        executor.session = session

        when:
        executor.createTaskMonitor()
        then:
        session.getExecConfigProp('array', 'target', 'local') >> 'nope'
        executorFactory.getExecutor('nope', session) >> new NopeExecutor()
        IllegalArgumentException e = thrown()
        e.message == "Executor 'nope' does not support array jobs"
    }

}
