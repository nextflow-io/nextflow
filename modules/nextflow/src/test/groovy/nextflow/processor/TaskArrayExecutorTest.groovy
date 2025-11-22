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

package nextflow.executor

import java.nio.file.Path

import nextflow.file.http.XPath
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskArrayExecutorTest extends Specification {

    static abstract class TestArrayExecutor implements TaskArrayExecutor {

    }

    @Unroll
    def 'should get array work dir' () {
        given:
        def task = Mock(TaskRun) {  getWorkDir() >> WORK_DIR }
        def handler = Mock(TaskHandler) { getTask()>>task }
        and:
        def executor = Spy(TestArrayExecutor) { isFusionEnabled() >> FUSION }
        when:
        def result = executor.getArrayWorkDir(handler)
        then:
        result == EXPECTED
        
        where:
        FUSION  | WORK_DIR                           |  EXPECTED
        false   | Path.of('/work/dir')               | '/work/dir'
        false   | XPath.get('http://foo.com/work')   | 'http://foo.com/work'
        true    | XPath.get('http://foo.com/work')   | '/fusion/http/foo.com/work'
    }

    @Unroll
    def 'should get array launch command' () {
        given:
        def executor = Spy(TestArrayExecutor) {
            isFusionEnabled() >> FUSION
            isWorkDirDefaultFS() >> !FUSION
        }

        expect:
        executor.getArrayLaunchCommand(WORK_DIR) == EXPECTED

        where:
        FUSION      | WORK_DIR              | EXPECTED
        false       | '/work/dir'           | 'bash /work/dir/.command.run 2>&1 > /work/dir/.command.log'
        true        | '/fusion/work/dir'    | 'bash /fusion/work/dir/.command.run'
    }

}
