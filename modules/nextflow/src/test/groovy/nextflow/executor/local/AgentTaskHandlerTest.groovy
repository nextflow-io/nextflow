/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.spi.FileSystemProvider

import nextflow.Session
import nextflow.file.FileHolder
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.ArrayBag
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The in-JVM half of agent path-input parity: a native agent task runs no wrapper script, so the
 * handler is what puts the declared input files into the work dir. Driven against a hand-built
 * {@link TaskRun} because the script-level harness substitutes every executor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentTaskHandlerTest extends Specification {

    @TempDir
    Path tempDir

    def 'should materialize the declared input files under their stage name'() {
        given:
        final task = taskWith(['reads.fq': sourceFile('reads.fq')])

        when:
        AgentTaskHandler.materializeInputs(task)

        then:
        final staged = task.workDir.resolve('reads.fq')
        Files.isSymbolicLink(staged)
        staged.text == 'DATA'
    }

    def 'should honour stageInMode'() {
        given:
        final source = sourceFile('reads.fq')
        final task = taskWith(['reads.fq': source], [stageInMode: mode])

        when:
        AgentTaskHandler.materializeInputs(task)

        then:
        final staged = task.workDir.resolve('reads.fq')
        Files.isSymbolicLink(staged) == symlink
        staged.text == 'DATA'
        and: 'a rellink resolves through a relative target, an absolute symlink through an absolute one'
        !symlink || Files.readSymbolicLink(staged).isAbsolute() == absoluteTarget

        where:
        mode      || symlink | absoluteTarget
        null      || true    | true
        'symlink' || true    | true
        'rellink' || true    | false
        'link'    || false   | true      // hard link: `absoluteTarget` is unused
        'copy'    || false   | true      // ditto
    }

    def 'should reject an unknown stage-in mode, exactly as a scriptlet task does'() {
        given:
        final task = taskWith(['reads.fq': sourceFile('reads.fq')], [stageInMode: 'teleport'])

        when:
        AgentTaskHandler.materializeInputs(task)

        then: 'the same message SimpleFileCopyStrategy raises -- a typo means one thing, not two'
        def e = thrown(IllegalArgumentException)
        e.message == 'Unknown stage-in strategy: teleport'
    }

    def 'should skip materialization when the work dir is not on the local filesystem'() {
        given: 'a work dir on a foreign provider, which cannot create links'
        final provider = Mock(FileSystemProvider) { getScheme() >> 's3' }
        final fs = Mock(FileSystem) { provider() >> provider }
        final workDir = Mock(Path) { getFileSystem() >> fs }
        and:
        final task = taskWith(['reads.fq': sourceFile('reads.fq')])
        task.workDir = workDir

        when:
        AgentTaskHandler.materializeInputs(task)

        then: 'no link attempt is made, so no UnsupportedOperationException escapes submit()'
        noExceptionThrown()
        0 * provider.createSymbolicLink(*_)
        0 * provider.createLink(*_)
        0 * workDir.resolve(_)
    }

    def 'should materialize a stage name carrying a sub-directory'() {
        given:
        final task = taskWith(['sub/dir/reads.fq': sourceFile('reads.fq')])

        when:
        AgentTaskHandler.materializeInputs(task)

        then:
        task.workDir.resolve('sub/dir/reads.fq').text == 'DATA'
    }

    def 'should replace a stale entry left by a previous attempt'() {
        given:
        final task = taskWith(['reads.fq': sourceFile('reads.fq')])
        and: 'a leftover regular file under the same name'
        task.workDir.resolve('reads.fq').text = 'STALE'

        when:
        AgentTaskHandler.materializeInputs(task)

        then:
        task.workDir.resolve('reads.fq').text == 'DATA'
    }

    def 'should do nothing when the task declares no input files'() {
        given:
        final task = taskWith([:])

        when:
        AgentTaskHandler.materializeInputs(task)

        then:
        noExceptionThrown()
        task.workDir.toFile().list().length == 0
    }

    def 'the agent executor dispatches a native task to the staging handler'() {
        given:
        final executor = new AgentExecutor(session: Mock(Session))
        final task = taskWith([:])
        task.type = ScriptType.GROOVY

        expect:
        executor.createTaskHandler(task) instanceof AgentTaskHandler
    }

    // -----------------------------------------------------------------------

    private Path sourceFile(String name) {
        final dir = Files.createTempDirectory(tempDir, 'src')
        final file = dir.resolve(name)
        file.text = 'DATA'
        return file
    }

    private TaskRun taskWith(Map<String,Path> inputs, Map directives = [:]) {
        final task = new TaskRun(id: null)
        task.workDir = Files.createTempDirectory(tempDir, 'work')
        task.config = new TaskConfig(directives)
        task.inputFiles = new ArrayBag<>(inputs.collect { name, path ->
            new FileHolder(path).withName(name)
        })
        return task
    }

}
