/*
 * Copyright 2020-2022, Seqera Labs
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

import java.nio.file.Paths

import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GridTaskHandlerTest extends Specification {

    def 'should check retry predicate' () {
        given:
        def handler = new GridTaskHandler()

        when:
        def predicate = handler.retryCondition("Socket timed out")
        then:
        predicate.test(new ProcessNonZeroExitStatusException('Error', 'Socket timed out', 1, []))
        and:
        predicate.test(new ProcessNonZeroExitStatusException('Error', 'error\nBatch job submission failed\nSocket timed out on send/recv operation', 1, [] ))
        and:
        !predicate.test(new ProcessNonZeroExitStatusException('Error', 'OK', 0, []))

    }

    def 'should capture error cause' () {
        given:
        def task = new TaskRun(name: 'foo', workDir: Paths.get('/some/work'))
        def exec = Mock(AbstractGridExecutor)
        def handler = Spy(new GridTaskHandler(task, exec))

        when:
        handler.submit()

        then:
        handler.safeExecute( _ )  >> { throw new ProcessNonZeroExitStatusException("Submit failed", "The limit is invalid", 10, ['qsub', 'foo']) }
        and:
        exec.createBashWrapperBuilder(task) >> Mock(BashWrapperBuilder)
        exec.pipeLauncherScript() >> false
        and:
        handler.createProcessBuilder() >> GroovyMock(ProcessBuilder)
        and:
        thrown(ProcessFailedException)
        and:
        task.stdout ==  "The limit is invalid"
        task.exitStatus == 10

    }
}
